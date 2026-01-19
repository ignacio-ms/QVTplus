function locEntry = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT)
    % Initialize parameters
    if length(correspondenceDict.(keyName)) == 1
        bestSegment = correspondenceDict.(keyName);
        segmentMask = data_struct.branchList(:, 4) == bestSegment;
        segmentIndices = find(segmentMask);
        segmentLength = length(segmentIndices);
        
        % Get MCA-PCA region for RCOMM/LCOMM proximity scoring
        mcaPcaRegion = [];
        if ismember(keyName, {'RCOMM', 'LCOMM'})
            mcaPcaRegion = getMCAPCARegion(keyName, correspondenceDict, data_struct);
        end
        
        % Use circularity-based selection if available
        if segmentLength >= 5 && isfield(data_struct, 'diam_val')
            % Consider points in the middle portion (avoid ends)
            startIdx = max(1, round(segmentLength * 0.25));
            endIdx = min(segmentLength, round(segmentLength * 0.75));
            candidateRange = startIdx:endIdx;
            
            bestScore = -inf;
            bestLOCIdx = [];
            
            for idx = candidateRange
                rowIdx = segmentIndices(idx);
                circularity = find_LOCs('getCircularity', data_struct, rowIdx);
                
                % For RCOMM/LCOMM, combine circularity with proximity to MCA-PCA region
                if ismember(keyName, {'RCOMM', 'LCOMM'}) && ~isempty(mcaPcaRegion)
                    candidatePos = data_struct.branchList(rowIdx, 1:3);
                    proximityScore = computeMCAPCAProximity(candidatePos, mcaPcaRegion);
                    % Combine: 60% circularity, 40% proximity
                    combinedScore = 0.25 * circularity + 0.75 * proximityScore;
                else
                    % For other vessels, use circularity only
                    combinedScore = circularity;
                end
                
                if combinedScore > bestScore
                    bestScore = combinedScore;
                    bestLOCIdx = idx;
                end
            end
            
            if ~isempty(bestLOCIdx)
                locEntry = [bestSegment, bestLOCIdx];
            else
                % Fallback
                if segmentLength == 5
                    LOC = 3;
                elseif segmentLength == 6
                    LOC = 4;
                else
                    LOC = 5;
                end
                locEntry = [bestSegment, LOC];
            end
        else
            % Fallback to original logic
            if segmentLength == 5
                LOC = 3;
            elseif segmentLength == 6
                LOC = 4;
            else
                LOC = 5;
            end
            locEntry = [bestSegment, LOC];
        end
        return;
    end

    segmentIndices = correspondenceDict.(keyName);
    bestSegment = [];
    highestZ = -inf;
    lowestZ = inf;
    largestY = -inf;
    lowestDist = inf;
    bestProximity = -inf; % For RCOMM/LCOMM MCA-PCA proximity scoring

    % Loop through each segment index and process it
    for segIdx = segmentIndices'
        segmentPositions = data_struct.branchList(data_struct.branchList(:, 4) == segIdx, 1:3);

        if size(segmentPositions, 1) > 4
            if ismember(keyName, {'RACA', 'LACA'})
                % Find segment with the highest Z value
                maxZ = mean(segmentPositions(:, 3));
                if maxZ > highestZ
                    bestSegment = segIdx;
                    highestZ = maxZ;
                end

            elseif ismember(keyName, {'RPCA', 'LPCA'})
                % % Prioritize based on Y and Z criteria
                % firstPoint = segmentPositions(1, :);
                % lastPoint = segmentPositions(end, :);
                % 
                % if lastPoint(2) + 5 < firstPoint(2) % Check Y condition
                %     meanY = mean(segmentPositions(:, 2));
                %     meanZ = mean(segmentPositions(:, 3));
                %     if meanY + 2 > largestY && meanZ - 15 < lowestZ
                %         bestSegment = segIdx;
                %         largestY = meanY;
                %         lowestZ = meanZ;
                %     end
                % end
                bestSegment = [];
                bestSignalQuality = -inf;
            
                % Merge PCA + PC2 candidates
                segmentSubIndices = correspondenceDict.(keyName);
                if strcmp(keyName, 'RPCA') && isfield(correspondenceDict, 'RPC2')
                    segmentSubIndices = [segmentSubIndices; correspondenceDict.('RPC2')];
                elseif strcmp(keyName, 'LPCA') && isfield(correspondenceDict, 'LPC2')
                    segmentSubIndices = [segmentSubIndices; correspondenceDict.('LPC2')];
                end

                segmentSubIndices = segmentSubIndices';
            
                for segIdx2 = segmentSubIndices
                    segmentMask = data_struct.branchList(:, 4) == segIdx2;
                    if nnz(segmentMask) < 5
                        continue;
                    end
            
                    flowVals = data_struct.flowPulsatile_val(segmentMask, :);  % [pixels x time]
                    if size(flowVals, 1) > 1
                        corrMatrix = corr(flowVals');  % transpose to [time x pixels], get [pixels x pixels] correlation
                        upperTri = triu(corrMatrix, 1);  % exclude diagonal and lower triangle
                        corrVals = upperTri(upperTri ~= 0);
                        signalStrength = mean(corrVals);  % or median(corrVals)
                    else
                        signalStrength = 0;  % Not enough data to compute correlation
                    end
            
                    if signalStrength > bestSignalQuality
                        bestSegment = segIdx2;
                        bestSignalQuality = signalStrength;
                    end
                end
            elseif ismember(keyName, {'RMCA', 'LMCA'})
                % Prioritize segments closer to centerline
                meanX = mean(segmentPositions(:, 1));
                centerlineDist = abs(meanX - (size(multiQVT, 1) / 2) - 3);
                if centerlineDist < lowestDist
                    bestSegment = segIdx;
                    lowestDist = centerlineDist;
                end
            elseif ismember(keyName, {'RCOMM', 'LCOMM'})
                % For RCOMM/LCOMM, prioritize segments within MCA-PCA region
                % Get MCA-PCA region boundaries
                mcaPcaRegion = getMCAPCARegion(keyName, correspondenceDict, data_struct);
                
                if ~isempty(mcaPcaRegion)
                    % Check if segment Y coordinate is within ACA-PCA Y range (in image space)
                    % In image space: Y increases top to bottom (anterior to posterior)
                    % Valid range: ACA_y_pos < y_pos < PCA_y_pos (in image coordinates)
                    % Region spans from minimum ACA Y to maximum PCA Y (all Y-axis stacks between them)
                    segmentMean = mean(segmentPositions, 1); % Mean position [X, Y, Z]
                    y = segmentMean(2);
                    withinYRegion = (y > mcaPcaRegion.yMin && y < mcaPcaRegion.yMax);
                    
                    if withinYRegion
                        % Only consider segments within the Y region
                        segmentProximity = computeMCAPCAProximity(segmentMean, mcaPcaRegion);
                        
                        % Use proximity score (higher = better, within region)
                        if isempty(bestSegment) || segmentProximity > bestProximity
                            bestSegment = segIdx;
                            bestProximity = segmentProximity;
                        end
                    end
                    % Discard segments outside the Y region (do nothing, continue to next segment)
                else
                    % Fallback: use Z position if MCA/PCA not available
                    meanZ = mean(segmentPositions(:, 3));
                    if meanZ > highestZ
                        bestSegment = segIdx;
                        highestZ = meanZ;
                    end
                end
            end
        end
    end

    % Check and update for RPC2/LPC2
    % if strcmp(keyName, 'RPCA') && ~isempty(bestSegment) && isfield(correspondenceDict, 'RPC2')
    %     bestSegment = checkSecondarySegments(bestSegment, correspondenceDict.('RPC2'), data_struct);
    % elseif strcmp(keyName, 'LPCA') && ~isempty(bestSegment) && isfield(correspondenceDict, 'LPC2')
    %     bestSegment = checkSecondarySegments(bestSegment, correspondenceDict.('LPC2'), data_struct);
    % end

    % Assign final segment or default to NaN
    if ~isempty(bestSegment)
        % Get all points in the selected segment
        segmentMask = data_struct.branchList(:, 4) == bestSegment;
        segmentIndices = find(segmentMask);
        segmentLength = length(segmentIndices);
        
            % Get MCA-PCA region for RCOMM/LCOMM proximity scoring
            mcaPcaRegion = [];
            if ismember(keyName, {'RCOMM', 'LCOMM'})
                mcaPcaRegion = getMCAPCARegion(keyName, correspondenceDict, data_struct);
            end
            
            % Evaluate circularity for candidate points along the segment
        if segmentLength >= 5 && isfield(data_struct, 'diam_val')
            % Consider points in the middle portion of the segment (avoid ends)
            % Evaluate middle 60% of segment (avoid first/last 20%)
            startIdx = max(1, round(segmentLength * 0.25));
            endIdx = min(segmentLength, round(segmentLength * 0.75));
            candidateRange = startIdx:endIdx;
            
            bestScore = -inf;
            bestLOCIdx = [];
            
            for idx = candidateRange
                rowIdx = segmentIndices(idx);
                candidatePos = data_struct.branchList(rowIdx, 1:3);
                
                % For RCOMM/LCOMM, discard candidates outside ACA-PCA Y region
                if ismember(keyName, {'RCOMM', 'LCOMM'}) && ~isempty(mcaPcaRegion)
                    % Check if candidate Y coordinate is within ACA-PCA Y range (in image space)
                    % In image space: Y increases top to bottom (anterior to posterior)
                    % Valid range: ACA_y_pos < y_pos < PCA_y_pos (in image coordinates)
                    % Region spans from minimum ACA Y to maximum PCA Y (all Y-axis stacks between them)
                    y = candidatePos(2);
                    withinYRegion = (y > mcaPcaRegion.yMin && y < mcaPcaRegion.yMax);
                    
                    if ~withinYRegion
                        % Discard candidates outside the Y region
                        continue;
                    end
                    
                    % Candidate is within Y region, compute scores
                    circularity = find_LOCs('getCircularity', data_struct, rowIdx);
                    proximityScore = computeMCAPCAProximity(candidatePos, mcaPcaRegion);
                    % Combine: 25% circularity, 75% proximity
                    combinedScore = 0.25 * circularity + 0.75 * proximityScore;
                else
                    % For other vessels, use circularity only
                    circularity = find_LOCs('getCircularity', data_struct, rowIdx);
                    combinedScore = circularity;
                end
                
                if combinedScore > bestScore
                    bestScore = combinedScore;
                    bestLOCIdx = idx;  % Position index within branch
                end
            end
            
            if ~isempty(bestLOCIdx)
                locEntry = [bestSegment, bestLOCIdx];
            else
                % For RCOMM/LCOMM, if no valid candidates within Y range, don't assign LOC
                if ismember(keyName, {'RCOMM', 'LCOMM'}) && ~isempty(mcaPcaRegion)
                    % No candidates found within valid Y range: PCA_y_pos < y_pos < ACA_y_pos
                    % Return NaN to indicate LOC should not exist
                    locEntry = [NaN, NaN];
                else
                    % Fallback to original logic for other vessels if no valid candidates
                    if segmentLength == 5
                        LOC = 3;
                    elseif segmentLength == 6
                        LOC = 4;
                    else
                        LOC = 5;
                    end
                    locEntry = [bestSegment, LOC];
                end
            end
        else
            % Very short segment or no circularity data, use original logic
            if segmentLength == 5
                LOC = 3;
            elseif segmentLength == 6
                LOC = 4;
            else
                LOC = 5;
            end
            locEntry = [bestSegment, LOC];
        end
    else
        locEntry = [NaN, NaN];
    end
end

function region = getMCAPCARegion(keyName, correspondenceDict, data_struct)
    % Get the Y-axis region between ACA and PCA vessels for RCOMM/LCOMM candidate filtering
    % Returns a structure with Y range: {yMin, yMax}
    % In image space: Y increases from top to bottom (anterior to posterior)
    % yMin = minimum Y coordinate of ACA (anterior, smaller Y in image space) - lower bound
    % yMax = maximum Y coordinate of PCA (posterior, larger Y in image space) - upper bound
    % Valid range in image coordinates: ACA_y_pos < y_pos < PCA_y_pos
    % Empty if PCA or ACA not available
    
    region = [];
    
    % Determine which ACA/PCA to use based on keyName
    if strcmp(keyName, 'RCOMM')
        acaKey = 'RACA';
        pcaKey = 'RPCA';
    elseif strcmp(keyName, 'LCOMM')
        acaKey = 'LACA';
        pcaKey = 'LPCA';
    else
        return; % Not a COMM vessel
    end
    
    % Get PCA Y positions (lower bound - minimum Y coordinate)
    pcaY = [];
    if isfield(correspondenceDict, pcaKey) && ~isempty(correspondenceDict.(pcaKey))
        pcaSegments = correspondenceDict.(pcaKey);
        for segID = pcaSegments(:)'
            pcaMask = data_struct.branchList(:, 4) == segID;
            if any(pcaMask)
                pcaPositions = data_struct.branchList(pcaMask, 1:3);
                pcaY = [pcaY; pcaPositions(:, 2)]; % Y coordinates
            end
        end
    end
    
    % Get ACA Y positions (upper bound - maximum Y coordinate)
    acaY = [];
    if isfield(correspondenceDict, acaKey) && ~isempty(correspondenceDict.(acaKey))
        acaSegments = correspondenceDict.(acaKey);
        for segID = acaSegments(:)'
            acaMask = data_struct.branchList(:, 4) == segID;
            if any(acaMask)
                acaPositions = data_struct.branchList(acaMask, 1:3);
                acaY = [acaY; acaPositions(:, 2)]; % Y coordinates
            end
        end
    end
    
    % If both PCA and ACA are available, define Y region between them
    % In image space: Y increases from top to bottom (anterior to posterior)
    % ACA is anterior (smaller Y), PCA is posterior (larger Y)
    % Valid range for RCOMM/LCOMM: ACA_y_pos < y_pos < PCA_y_pos (in image coordinates)
    if ~isempty(pcaY) && ~isempty(acaY)
        % Y region: from minimum ACA Y to maximum PCA Y (in image space)
        % This covers all "stacks of the y axis" between ACA and PCA
        yMin = min(acaY); % Minimum ACA Y position (anterior, smaller Y in image space)
        yMax = max(pcaY); % Maximum PCA Y position (posterior, larger Y in image space)
        
        region = struct('yMin', yMin, 'yMax', yMax);
    end
end

function proximityScore = computeMCAPCAProximity(position, region)
    % Compute proximity score for a position relative to PCA-ACA Y region
    % Returns 1.0 if within Y region, decreasing to 0.0 as distance increases
    % position: [X, Y, Z] coordinates
    % region: struct with {yMin, yMax} - Y-axis range from minimum PCA Y to maximum ACA Y
    % Valid range: PCA_y_pos < y_pos < ACA_y_pos
    
    if isempty(region) || numel(position) < 3
        proximityScore = 0.0;
        return;
    end
    
    y = position(2);
    
    % Check if within Y region (should already be filtered, but compute score for ranking)
    % Valid range: PCA_y_pos < y_pos < ACA_y_pos
    if y > region.yMin && y < region.yMax
        % Within Y region: score = 1.0
        % Candidates within the Y range get full proximity score
        proximityScore = 1.0;
    else
        % Outside Y region: compute distance to Y boundaries (shouldn't happen after filtering)
        if y < region.yMin
            yDist = region.yMin - y;
        elseif y > region.yMax
            yDist = y - region.yMax;
        else
            yDist = 0;
        end
        
        % Normalize distance (use Y range as reference)
        yRange = region.yMax - region.yMin;
        
        if yRange > 0
            normalizedDist = yDist / yRange;
            % Score decreases exponentially with distance
            % Score = 1.0 at dist=0, ~0.5 at dist=yRange, ~0.0 at dist=2*yRange
            proximityScore = exp(-normalizedDist * 2);
        else
            proximityScore = 0.0;
        end
    end
end

function bestSegment = checkSecondarySegments(bestSegment, secondaryIndices, data_struct)
    largestY2 = -inf;

    for segIdx = secondaryIndices'
        segmentPositions = data_struct.branchList(data_struct.branchList(:, 4) == segIdx, 1:3);

        if size(segmentPositions, 1) < 5
            continue;
        end

        firstPoint = segmentPositions(1, :);
        lastPoint = segmentPositions(end, :);

        if lastPoint(2) + 5 < firstPoint(2) && (lastPoint(3) + 3) > firstPoint(3) % Conditions on Y and Z
            meanY = mean(segmentPositions(:, 2));
            if meanY > largestY2
                bestSegment = segIdx;
                largestY2 = meanY;
            end
        end
    end
end

