function varargout = find_LOCs(funcName, varargin)
    % Dispatcher function
    switch funcName
        case 'extractBranchInfo'
            varargout{1} = extractBranchInfo(varargin{:});
        case 'findMaxSingleZ'
            varargout{1} = findMaxSingleZ(varargin{:});
        case 'extractLocation'
            varargout{1} = extractLocation(varargin{:});
        case 'ensureMinZOffset'
            varargout{1} = ensureMinZOffset(varargin{:});
        case 'extractSSSV'
            varargout{1} = extractSSSV(varargin{:});
        case 'extractLTSV'
            varargout{1} = extractLTSV(varargin{:});
        case 'extractRTSV'
            varargout{1} = extractRTSV(varargin{:});
        case 'extractSTRV'
            varargout{1} = extractSTRV(varargin{:});
        case 'getCircularity'
            varargout{1} = getCircularity(varargin{:});
        case 'selectLOCByCircularity'
            varargout{1} = selectLOCByCircularity(varargin{:});
        otherwise
            error('Unknown function name: %s', funcName);
    end
end

function branchInfo = extractBranchInfo(data_struct, branchLabels)
    branchInfo = [];
    for i = 1:length(branchLabels)
        branchInfo = [branchInfo; data_struct.branchList(data_struct.branchList(:, 4) == branchLabels(i), :)];
    end
end

function maxZ = findMaxSingleZ(z_LICA, z_RICA, z_BA, info_LICA, info_RICA, info_BA)
    maxZ = -inf;
    for z = max([z_LICA; z_RICA; z_BA]):-1:min([z_LICA; z_RICA; z_BA])
        count_LICA = sum(info_LICA(:, 3) == z);
        count_RICA = sum(info_RICA(:, 3) == z);
        count_BA = sum(info_BA(:, 3) == z);

        if count_LICA == 1 && count_RICA == 1 && count_BA == 1
            maxZ = z;
            break;
        end
    end
end

function loc = extractLocation(info, zValue, zColumn)
    loc = info(find(info(:, zColumn) == zValue, 1), :);
end

function loc = ensureMinZOffset(info, loc, minOffset, topOffset)
    % Ensure minimum Z offset is at least minOffset
    if loc(1, 5) < minOffset
        loc(1, 5) = minOffset;
    end

    % Calculate max Z value for the same branch
    maxZ = max(info(info(:, 4) == loc(4), 5));

    % Adjust if the difference is less than minOffset
    if (maxZ - loc(1, 5)) < topOffset
        loc(1, 5) = maxZ - topOffset;
    end
end

function loc = extractSSSV(segment_ids, data_struct)
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments
    x = branchList(:,1);
    y = branchList(:,2);
    z = branchList(:,3);
    segment_id = branchList(:,4);

    dims = size(data_struct.segment);
    minX = min(x); maxX = max(x); 
    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);
    min_thrX = minX + 0.33 * (maxX - minX);
    max_thrX = maxX - 0.33 * (maxX - minX);
    posterior_thresh = minY + 0.4 * (maxY - minY);
    superior_thresh = minZ + 0.33 * (maxZ - minZ);
    % superior_mask  = z < midZ;
    posterior_mask = y <= posterior_thresh;
    superior_mask  = z >= superior_thresh;

    ss_mask = posterior_mask & superior_mask;
    subset = branchList(ss_mask, :);

    unique_segments = unique(subset(:,4));
    num_segments = length(unique_segments);

    if num_segments == 0
        loc = [];
        return;
    elseif num_segments <= 2
        counts = arrayfun(@(s) sum(subset(:,4) == s), unique_segments);
        [~, idx] = max(counts);
        ss_segment_id = unique_segments(idx);
    else
        lengths = zeros(num_segments,1);
        for i = 1:num_segments
            seg_id = unique_segments(i);
            lengths(i) = sum(subset(:,4) == seg_id);
        end
        [~, idx] = max(lengths);
        ss_segment_id = unique_segments(idx);
    end

    seg_points = branchList(segment_id == ss_segment_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);  % Return full row
    end

function loc = extractLTSV(segment_ids, data_struct)
    loc = extractLateralTSV(segment_ids, data_struct, 'left');
end

function loc = extractRTSV(segment_ids, data_struct)
    loc = extractLateralTSV(segment_ids, data_struct, 'right');
end

function loc = extractLateralTSV(segment_ids, data_struct, side)
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments

    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); val = branchList(:,4);

    minX = min(x); maxX = max(x);
    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);

    % x_left_thresh  = minX + 0.33 * (maxX - minX);
    % x_right_thresh = maxX - 0.33 * (maxX - minX);
    % posterior_thresh = minY + 0.33 * (maxY - minY);
    % z_inferior_thresh = minZ + 0.7 * (maxZ - minZ);
    x_left_thresh  = minX + 0.40 * (maxX - minX);
    x_right_thresh = maxX - 0.40 * (maxX - minX);
    posterior_thresh = minY + 0.4 * (maxY - minY);
    z_inferior_thresh = minZ + 0.4 * (maxZ - minZ);

    % inferior_mask = z >= z_inferior_thresh;
    inferior_mask = z <= z_inferior_thresh;
    posterior_mask = y <= posterior_thresh;

    if strcmp(side, 'left')
        % region_mask = x <= x_left_thresh & inferior_mask & posterior_mask;
        region_mask = x <= x_left_thresh & inferior_mask & posterior_mask;
    else
        % region_mask = x >= x_right_thresh & inferior_mask & posterior_mask;
        region_mask = x >= x_right_thresh & inferior_mask & posterior_mask;
    end

    subset = branchList(region_mask, :);
    segments = unique(subset(:,4));
    best_id = NaN;
    best_length = 0;
    z_std_thresh = 3.5;

    for i = 1:length(segments)
        seg_id = segments(i);
        seg_points = subset(subset(:,4) == seg_id, :);
        num_pts = size(seg_points,1);

        if num_pts <= 40
            z_std = std(seg_points(:,3));
            if z_std < z_std_thresh && num_pts > best_length
                best_length = num_pts;
                best_id = seg_id;
            end
        else
            quarter = floor(num_pts / 4);
            parts = {seg_points(1:quarter,:); seg_points(quarter+1:2*quarter,:); ...
                     seg_points(2*quarter+1:3*quarter,:); seg_points(3*quarter+1:end,:)};

            for part = parts
                p = part{1};
                if strcmp(side, 'left')
                    x_ok = p(:,1) <= x_left_thresh;
                else
                    x_ok = p(:,1) >= x_right_thresh;
                end
                y_ok = p(:,2) <= posterior_thresh;
                % z_ok = p(:,3) >= z_inferior_thresh;
                z_ok = p(:,3) <= z_inferior_thresh;
                region_mask = x_ok & y_ok & z_ok;

                if mean(region_mask) >= 0.9
                    z_std = std(p(:,3));
                    len = size(p,1);
                    if z_std < z_std_thresh && len > best_length
                        best_length = len;
                        best_id = seg_id;
                    end
                end
            end
        end
    end
    segment_id = branchList(:,4);  % Define it here
    seg_points = branchList(segment_id == best_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);
end

function loc = extractSTRV(segment_ids, data_struct, sssv_segment_id)
    % Extract STRV LOC, excluding segments already assigned to SSSV
    % 
    % Inputs:
    %   segment_ids - Segment IDs to consider for STRV
    %   data_struct - Data structure containing branchList
    %   sssv_segment_id (optional) - Segment ID already assigned to SSSV (if any)
    %
    % Output:
    %   loc - Location row from branchList for STRV, or empty if no valid segment found
    
    if nargin < 3
        sssv_segment_id = [];
    end
    
    branchList = data_struct.branchList;
    branchList = branchList(ismember(branchList(:,4), segment_ids), :);  % Restrict to relevant segments
    x = branchList(:,1); y = branchList(:,2); z = branchList(:,3); segment_id = branchList(:,4);

    minY = min(y); maxY = max(y);
    minZ = min(z); maxZ = max(z);
    posterior_thresh = minY + 0.5 * (maxY - minY);
    superior_thresh  = minZ + 0.33 * (maxZ - minZ);

    posterior_mask = y <= posterior_thresh;
    % superior_mask  = z <= superior_thresh;
    superior_mask  = z >= superior_thresh;
    straight_mask = posterior_mask & superior_mask;

    straight_subset = branchList(straight_mask, :);
    straight_segments = unique(straight_subset(:,4));

    % Exclude SSSV segment if it was already assigned
    if ~isempty(sssv_segment_id) && ismember(sssv_segment_id, straight_segments)
        straight_segments = straight_segments(straight_segments ~= sssv_segment_id);
    end

    direction_threshold = 0.85;
    min_points = 15;
    % expected = [0; 1; -1]; expected = expected / norm(expected);
    expected = [0; 1; 1]; expected = expected / norm(expected);
    best_score = -Inf;
    straight_sinus_id = NaN;

    for i = 1:length(straight_segments)
        seg_id = straight_segments(i);
        if seg_id == 1, continue; end
        seg_points = straight_subset(straight_subset(:,4) == seg_id, 1:3);
        if size(seg_points,1) < min_points, continue; end

        centered = seg_points - mean(seg_points);
        [~, ~, V] = svd(centered, 'econ');
        direction = V(:,1);
        alignment = abs(dot(direction, expected));

        if alignment > direction_threshold && alignment > best_score
            best_score = alignment;
            straight_sinus_id = seg_id;
        end
    end

    if isnan(straight_sinus_id)
        loc = [];
        return;
    end

    seg_points = branchList(segment_id == straight_sinus_id, :);
    center = mean(seg_points(:,1:3), 1);
    [~, idx] = min(sum((seg_points(:,1:3) - center).^2, 2));
    loc = seg_points(idx, :);
end

function circularity = getCircularity(data_struct, rowIdx)
    % Get circularity value for a given branchList row index
    % Returns diam_val (Rin^2/Rout^2) which is already computed
    % Value of 1 = perfect circle, <1 = irregular shape
    %
    % Inputs:
    %   data_struct - Data structure containing diam_val
    %   rowIdx - Row index in branchList (1-based)
    %
    % Output:
    %   circularity - Circularity metric (0 if invalid)
    
    if nargin < 2 || isempty(data_struct) || ~isfield(data_struct, 'diam_val')
        circularity = 0;
        return;
    end
    
    if rowIdx > 0 && rowIdx <= length(data_struct.diam_val)
        circularity = data_struct.diam_val(rowIdx);
        % Handle invalid values
        if ~isfinite(circularity) || circularity <= 0
            circularity = 0;
        end
    else
        circularity = 0;
    end
end

function bestLOC = selectLOCByCircularity(candidates, data_struct, fullInfo, minOffset, topOffset)
    % Select best LOC from candidates based on circularity and flow, avoiding extremes
    % 
    % Inputs:
    %   candidates - Candidate points (rows from branchList, N x M matrix)
    %   data_struct - Data structure with diam_val, flowPerHeartCycle_val, and branchList
    %   fullInfo - Full branch info for ensureMinZOffset
    %   minOffset, topOffset - Parameters for ensureMinZOffset
    %
    % Output:
    %   bestLOC - Best candidate point (full row) or empty if no valid candidates
    
    if isempty(candidates)
        bestLOC = [];
        return;
    end
    
    % Parameters for avoiding extremes (exclude first/last 20% of centerline)
    extremeThreshold = 0.2;
    
    % Pre-compute segment information to avoid extremes
    segmentInfo = struct();
    for i = 1:size(candidates, 1)
        candidate = candidates(i, :);
        if size(candidate, 2) < 5
            continue;
        end
        segID = candidate(4);
        if ~isfield(segmentInfo, sprintf('seg_%d', segID))
            % Get all points in this segment
            segMask = data_struct.branchList(:, 4) == segID;
            segClIndices = data_struct.branchList(segMask, 5);
            if ~isempty(segClIndices)
                minClIdx = min(segClIndices);
                maxClIdx = max(segClIndices);
                segLength = maxClIdx - minClIdx + 1;
                extremeRange = round(segLength * extremeThreshold);
                segmentInfo.(sprintf('seg_%d', segID)) = struct(...
                    'minClIdx', minClIdx, ...
                    'maxClIdx', maxClIdx, ...
                    'minValidClIdx', minClIdx + extremeRange, ...
                    'maxValidClIdx', maxClIdx - extremeRange);
            end
        end
    end
    
    % Get flow values for normalization
    validCandidates = [];
    flowValues = [];
    circularityValues = [];
    isExtremeFlags = [];
    
    for i = 1:size(candidates, 1)
        candidate = candidates(i, :);
        
        if size(candidate, 2) < 5
            continue; % Invalid candidate row
        end
        
        % Apply ensureMinZOffset constraint
        candidate = find_LOCs('ensureMinZOffset', fullInfo, candidate, minOffset, topOffset);
        
        % Find row index in branchList
        segID = candidate(4);
        clIdx = candidate(5);
        rowIdx = find(data_struct.branchList(:,4) == segID & ...
                     data_struct.branchList(:,5) == clIdx, 1);
        
        if isempty(rowIdx)
            continue;
        end
        
        % Check if candidate is at extreme of centerline
        isExtreme = false;
        segKey = sprintf('seg_%d', segID);
        if isfield(segmentInfo, segKey)
            segInfo = segmentInfo.(segKey);
            if clIdx < segInfo.minValidClIdx || clIdx > segInfo.maxValidClIdx
                isExtreme = true;
            end
        end
        
        % Get circularity
        circularity = getCircularity(data_struct, rowIdx);
        
        % Get flow value
        flow = 0;
        if isfield(data_struct, 'flowPerHeartCycle_val') && ...
           rowIdx > 0 && rowIdx <= length(data_struct.flowPerHeartCycle_val)
            flow = data_struct.flowPerHeartCycle_val(rowIdx);
            if ~isfinite(flow) || flow < 0
                flow = 0;
            end
        end
        
        % Store candidate information
        validCandidates = [validCandidates; candidate];
        circularityValues = [circularityValues; circularity];
        flowValues = [flowValues; flow];
        isExtremeFlags = [isExtremeFlags; isExtreme];
    end
    
    if isempty(validCandidates)
        % Fallback to first candidate
        if ~isempty(candidates)
            bestLOC = candidates(1, :);
            if size(bestLOC, 2) >= 5
                bestLOC = find_LOCs('ensureMinZOffset', fullInfo, bestLOC, minOffset, topOffset);
            end
        else
            bestLOC = [];
        end
        return;
    end
    
    % Normalize circularity and flow values for scoring
    % Circularity: higher is better (already normalized 0-1 typically)
    validCircularity = circularityValues(isfinite(circularityValues) & circularityValues > 0);
    if ~isempty(validCircularity)
        maxCircularity = max(validCircularity);
        if maxCircularity > 0
            normalizedCircularity = circularityValues / maxCircularity;
        else
            normalizedCircularity = zeros(size(circularityValues));
        end
    else
        normalizedCircularity = zeros(size(circularityValues));
    end
    
    % Flow: higher is better, normalize to 0-1
    validFlow = flowValues(isfinite(flowValues) & flowValues > 0);
    if ~isempty(validFlow)
        maxFlow = max(validFlow);
        if maxFlow > 0
            normalizedFlow = flowValues / maxFlow;
        else
            normalizedFlow = zeros(size(flowValues));
        end
    else
        normalizedFlow = zeros(size(flowValues));
    end
    
    % Combined score: weighted combination of circularity and flow
    % Penalize extreme positions
    % Weights: 60% circularity, 40% flow (adjustable)
    circularityWeight = 0.6;
    flowWeight = 0.4;
    
    combinedScore = circularityWeight * normalizedCircularity + flowWeight * normalizedFlow;
    
    % Apply penalty for extreme positions (reduce score by 50%)
    extremePenalty = 0.5;
    % combinedScore(isExtremeFlags) = combinedScore(isExtremeFlags) * extremePenalty;
    for i = 1:length(isExtremeFlags)
        if isExtremeFlags(i)
            combinedScore(i) = combinedScore(i) * extremePenalty;
        end
    end
    
    % Find best candidate
    [~, bestIdx] = max(combinedScore);
    bestLOC = validCandidates(bestIdx, :);
    
    % Final fallback if still no valid candidate
    if isempty(bestLOC) && ~isempty(candidates)
        bestLOC = candidates(1, :);
        if size(bestLOC, 2) >= 5
            bestLOC = find_LOCs('ensureMinZOffset', fullInfo, bestLOC, minOffset, topOffset);
        end
    end
end
