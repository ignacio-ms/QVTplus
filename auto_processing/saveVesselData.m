function saveVesselData(LOCs, data_struct, output_path)
    fields = fieldnames(LOCs);
    
    generateSummaryCenterline(LOCs, data_struct, output_path)
    for i = 1:numel(fields)
        vesselName = fields{i};

        % RCOMM and LCOMM are now handled separately (split from COMM based on RAS orientation)
        % Skip old 'COMM' if it exists (should not exist after splitting, but check for compatibility)
        if strcmp(vesselName, 'COMM')
            continue;
        end

        vesselInfo = LOCs.(vesselName);
        
        % Check if vessel info is valid (not NaN, not empty)
        if isempty(vesselInfo) || numel(vesselInfo) < 2 || ...
           isnan(vesselInfo(1)) || isnan(vesselInfo(2)) || ...
           ~isfinite(vesselInfo(1)) || ~isfinite(vesselInfo(2))
            % Skip invalid vessels (e.g., missing RCOMM/LCOMM)
            continue;
        end
        
        vesselNumber = vesselInfo(1);
        pointOfInterest = vesselInfo(2);

        % Extract data and save
        saveVesselExcelData(vesselName, vesselNumber, pointOfInterest, data_struct, output_path);
        saveVesselImages(vesselName, vesselNumber, pointOfInterest, data_struct, output_path);
    end
end

function generateSummaryCenterline(LOCs, data_struct, output_path)
    % Vessel Labels
    vesselLabels = {
        'Left ICA', 'LICA';
        'Right ICA', 'RICA';
        'Basilar', 'BASI';
        'Left MCA', 'LMCA';
        'Right MCA', 'RMCA';
        'Left PCA', 'LPCA';
        'Right PCA', 'RPCA';
        'Left ACA', 'LACA';
        'Right ACA', 'RACA';
        'Right Communicating', 'RCOMM';
        'Left Communicating', 'LCOMM';
        'Straight Sinus', 'STRV';
        'Sagital Sinus', 'SSSV';
        'Left Transverse', 'LTSV';
        'Right Transverse', 'RTSV'
    };

    nBranchList = size(data_struct.branchList, 1);
    nFlow = numel(data_struct.flowPerHeartCycle_val);
    nPI = numel(data_struct.PI_val);
    fprintf('[generateSummaryCenterline] branchList rows = %d, flowPerHeartCycle_val length = %d, PI_val length = %d\n', ...
        nBranchList, nFlow, nPI);

    % Initialize the table
    summaryData = cell(size(vesselLabels, 1), 7);
    columnNames = {'Vessel Label', 'Centerline', 'Notes', 'Max Vel < 700 cm/s', ...
                   'Mean Flow ml/s', 'Pulsatility Index', 'Branch Number'};

    % Fill the table
    for i = 1:size(vesselLabels, 1)
        vesselLabel = vesselLabels{i, 1};
        locKey = vesselLabels{i, 2};
        summaryData{i, 1} = vesselLabel;

        if isfield(LOCs, locKey)
            vesselLOC = LOCs.(locKey);
            if isempty(vesselLOC) || numel(vesselLOC) < 2
                fprintf('[generateSummaryCenterline] %s: LOC has < 2 elements -> NaN\n', locKey);
                summaryData(i, :) = {vesselLabel, NaN, NaN, NaN, NaN, NaN, NaN};
            else
            summaryData{i, 2} = vesselLOC(2); % Centerline
            summaryData{i, 7} = vesselLOC(1); % Branch Number
            % LOC is (segment_id, centerline_point_index) from ICA/BA/main vessels, or (segment_id, row_index) from resolveLongVenousSegment
            rowIndex = find(data_struct.branchList(:, 4) == vesselLOC(1) & data_struct.branchList(:, 5) == vesselLOC(2), 1);
            if isempty(rowIndex) && vesselLOC(2) >= 1 && vesselLOC(2) <= nBranchList
                % Fallback: vesselLOC(2) may be stored as row index (e.g. venous from resolveLongVenousSegment)
                if data_struct.branchList(vesselLOC(2), 4) == vesselLOC(1)
                    rowIndex = vesselLOC(2);
                    fprintf('[generateSummaryCenterline] %s: used row index fallback (seg=%d, stored value=%d)\n', ...
                        locKey, vesselLOC(1), vesselLOC(2));
                end
            end
            if ~isempty(rowIndex)
                if rowIndex > nFlow || rowIndex > nPI
                    fprintf('[generateSummaryCenterline] %s: rowIndex=%d out of range (flow len=%d, PI len=%d) -> NaN\n', ...
                        locKey, rowIndex, nFlow, nPI);
                    summaryData{i, 4} = NaN;
                    summaryData{i, 5} = NaN;
                    summaryData{i, 6} = NaN;
                else
                    summaryData{i, 4} = data_struct.maxVel_val(rowIndex) < 700;
                    summaryData{i, 5} = data_struct.flowPerHeartCycle_val(rowIndex);
                    summaryData{i, 6} = data_struct.PI_val(rowIndex);
                end
            else
                fprintf('[generateSummaryCenterline] %s: no row match for (seg=%d, centerline/pt=%d); branchList seg range for seg: %s\n', ...
                    locKey, vesselLOC(1), vesselLOC(2), mat2str(unique(data_struct.branchList(data_struct.branchList(:,4)==vesselLOC(1), 5))));
                summaryData{i, 4} = NaN;
                summaryData{i, 5} = NaN;
                summaryData{i, 6} = NaN;
            end
            end
        else
            fprintf('[generateSummaryCenterline] %s: no LOC in struct -> all NaN\n', locKey);
            summaryData(i, :) = {vesselLabel, NaN, NaN, NaN, NaN, NaN, NaN};
        end
    end

    % Write to Excel using xlwrite (consistent with rest of codebase)
    % First write column headers
    xlwrite(fullfile(output_path, 'SummaryParamTool.xls'), columnNames, 'Summary_Centerline', 'A1');
    % Then write data starting from row 2
    xlwrite(fullfile(output_path, 'SummaryParamTool.xls'), summaryData, 'Summary_Centerline', 'A2');
end

function saveVesselExcelData(vesselName, vesselNumber, pointOfInterest, data_struct, output_path)
    % Extract relevant data
    branchList = data_struct.branchList;
    Logical_branch = branchList(:, 4) == vesselNumber;
    indicesInBranch = find(Logical_branch);
    
    % pointOfInterest can be either:
    % 1. A position index (1-based) within the branch (from processMainVessels)
    % 2. A row index in the full branchList (from resolveLongVenousSegment)
    % Check if pointOfInterest is a valid row index in the branch
    if pointOfInterest <= length(branchList) && Logical_branch(pointOfInterest)
        % It's a row index in the full branchList - convert to position index
        selectedPointIndex = pointOfInterest;
        pointPosition = find(indicesInBranch == pointOfInterest, 1);
        if isempty(pointPosition)
            pointPosition = ceil(length(indicesInBranch) / 2);
            selectedPointIndex = indicesInBranch(pointPosition);
        end
    elseif pointOfInterest > 0 && pointOfInterest <= length(indicesInBranch)
        % It's a position index within the branch
        selectedPointIndex = indicesInBranch(pointOfInterest);
        pointPosition = pointOfInterest;
    else
        % Invalid index - use the middle point as fallback
        warning('Point of interest %d is invalid for branch %d (vessel %s). Using middle point instead.', ...
                pointOfInterest, vesselNumber, vesselName);
        pointPosition = ceil(length(indicesInBranch) / 2);
        selectedPointIndex = indicesInBranch(pointPosition);
    end

    % Define a 5-point range centered on the selected point
    index_range = max(selectedPointIndex - 2, 1):min(selectedPointIndex + 2, length(branchList));
    index_range(~Logical_branch(index_range)) = [];
    
    % Get actual number of points in the range
    numPoints = length(index_range);

    % Time-averaged data calculations
    area = data_struct.area_val(index_range);
    area = [area; mean(area); std(area)];
    diam = data_struct.diam_val(index_range);
    diam = [diam; mean(diam); std(diam)];
    flowPerHeartCycle = data_struct.flowPerHeartCycle_val(index_range);
    flowPerHeartCycle = [flowPerHeartCycle; mean(flowPerHeartCycle); std(flowPerHeartCycle)];
    PI = data_struct.PI_val(index_range);
    PI = [PI; mean(PI); std(PI)];
    maxVel = data_struct.maxVel_val(index_range);
    maxVel = [maxVel; mean(maxVel); std(maxVel)];
    meanVel = data_struct.velMean_val(index_range);
    meanVel = [meanVel; mean(meanVel); std(meanVel)];
    RI = data_struct.RI_val(index_range);
    RI = [RI; mean(RI); std(RI)];

    % Time-resolved data
    flowPulsatile = data_struct.flowPulsatile_val(index_range, :);
    flowPulsatile = [flowPulsatile; mean(flowPulsatile, 1); std(flowPulsatile, 1)];

    % Labels for time-averaged data (use position index along vessel)
    % Create Labels to match the actual number of points in index_range
    % Find the position indices for each point in index_range
    Labels = zeros(1, numPoints);
    for i = 1:numPoints
        idx = index_range(i);
        pos = find(indicesInBranch == idx, 1);
        if ~isempty(pos)
            Labels(i) = pos;
        else
            Labels(i) = 0; % Invalid position
        end
    end
    Labels = [Labels, 0, 0]; % Adding placeholders for mean and std rows

    % Prepare Excel data for time-averaged data
    col_header = ({'Point along Vessel', 'Area (cm^2)', 'Area Ratio', 'Max Velocity (cm/s)',...
        'Mean Velocity (cm/s)','Average Flow(mL/s)','Pulsatility Index','Resistivity Index'});
    time_avg = vertcat(col_header,num2cell(real(horzcat(Labels',...
        area,diam,maxVel,meanVel,flowPerHeartCycle,PI,RI))));
    time_avg{end-1,1} = 'Mean';
    time_avg{end,1} = 'Standard Deviation';
    xlwrite([output_path filesep 'SummaryParamTool.xls'],time_avg,[vesselName '_T_averaged']);
    
    % save time-resolved
    spaces = repmat({''},1,data_struct.nframes-1);
    col_header2 = ({'Cardiac Time (ms)'});
    col_header3 = horzcat({'Point along Vessel','Flow (mL/s)'},spaces);
    col_header2 = horzcat(col_header2, num2cell(real(data_struct.timeres/1000*linspace(1,data_struct.nframes,data_struct.nframes))));
    time_resolve = vertcat(col_header2, col_header3, num2cell(real(horzcat(Labels',flowPulsatile))));
    time_resolve{end-1,1} = 'Mean';
    time_resolve{end,1} = 'Standard Deviation';
    xlwrite([output_path filesep 'SummaryParamTool.xls'],time_resolve,[vesselName '_T_resolved']);
end


function saveVesselImages(vesselName, vesselNumber, pointOfInterest, data_struct, output_path)
    % Extract relevant data
    branchList = data_struct.branchList;
    Logical_branch = branchList(:, 4) == vesselNumber;
    indicesInBranch = find(Logical_branch);
    
    % pointOfInterest can be either:
    % 1. A position index (1-based) within the branch (from processMainVessels)
    % 2. A row index in the full branchList (from resolveLongVenousSegment)
    % Check if pointOfInterest is a valid row index in the branch
    if pointOfInterest <= length(branchList) && Logical_branch(pointOfInterest)
        % It's a row index in the full branchList
        selectedPointIndex = pointOfInterest;
    elseif pointOfInterest > 0 && pointOfInterest <= length(indicesInBranch)
        % It's a position index within the branch
        selectedPointIndex = indicesInBranch(pointOfInterest);
    else
        % Invalid index - use the middle point as fallback
        warning('Point of interest %d is invalid for branch %d (vessel %s). Using middle point instead.', ...
                pointOfInterest, vesselNumber, vesselName);
        midPoint = ceil(length(indicesInBranch) / 2);
        selectedPointIndex = indicesInBranch(midPoint);
    end

    % Define a 5-point range centered on the selected point
    index_range = max(selectedPointIndex - 2, 1):min(selectedPointIndex + 2, length(branchList));
    index_range(~Logical_branch(index_range)) = [];

    % Cross-sectional image processing
    imdim = sqrt(size(data_struct.segmentFull, 2));
    BranchSlice = data_struct.segmentFull(index_range, :);
    cdSlice = data_struct.timeMIPcrossection(index_range, :);
    velSlice = data_struct.vTimeFrameave(index_range, :);
    magSlice = data_struct.MAGcrossection(index_range, :);

    subL = size(BranchSlice, 1);
    FinalImage = zeros(imdim, imdim, 1, 4 * subL);
    temp = 1;

    % Generate montage for cross-sections
    for q = 1:subL
        CDcross = reshape(cdSlice(q, :), imdim, imdim) ./ max(cdSlice(q, :));
        Vcross = reshape(velSlice(q, :), imdim, imdim) ./ max(velSlice(q, :));
        Magcross = reshape(magSlice(q, :), imdim, imdim) ./ max(magSlice(q, :));
        Maskcross = reshape(BranchSlice(q, :), imdim, imdim);

        FinalImage(:, :, 1, temp) = Magcross;
        FinalImage(:, :, 1, temp + 1) = CDcross;
        FinalImage(:, :, 1, temp + 2) = Vcross;
        FinalImage(:, :, 1, temp + 3) = Maskcross;
        temp = temp + 4;
    end

    % Save montage as an image
    montage(FinalImage, 'Size', [subL, 4]);
    savename = sprintf('%s_Vessel_%d_Point_%d_Slicesview.jpg', vesselName, vesselNumber, pointOfInterest);
    saveas(gcf, fullfile(output_path, savename));
    close(gcf);
end
