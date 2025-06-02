function saveVesselData(LOCs, data_struct, output_path)
    fields = fieldnames(LOCs);
    
    generateSummaryCenterline(LOCs, data_struct, output_path)
    for i = 1:numel(fields)
        vesselName = fields{i};

        % TODO: estamos saltando la COMM porque habria que tratarla distinto para
        % seleccionar el punto de interes (o seleccionar todos sus puntos para la media)
        if strcmp(vesselName, 'COMM')
            continue;
        end

        vesselInfo = LOCs.(vesselName);
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
        'Straight Sinus', 'STRV';
        'Sagital Sinus', 'SSSV';
        'Left Transverse', 'LTSV';
        'Right Transverse', 'RTSV'
    };

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
            % Extract LOC information
            vesselLOC = LOCs.(locKey);
            summaryData{i, 2} = vesselLOC(2); % Centerline
            summaryData{i, 7} = vesselLOC(1); % Branch Number
            % Compute other metrics
            rowIndex = find(data_struct.branchList(:, 4) == vesselLOC(1) & data_struct.branchList(:, 5) == vesselLOC(2));
            summaryData{i, 4} = data_struct.maxVel_val(rowIndex) < 700; % Max Vel < 700
            summaryData{i, 5} = data_struct.flowPerHeartCycle_val(rowIndex); % Mean Flow
            summaryData{i, 6} = data_struct.PI_val(rowIndex); % Pulsatility Index
        else
            summaryData(i, :) = {vesselLabel, NaN, NaN, NaN, NaN, NaN, NaN};
        end
    end

    % Write to Excel
    summaryTable = cell2table(summaryData, 'VariableNames', columnNames);
    writetable(summaryTable, fullfile(output_path, 'SummaryParamTool.xls'), 'Sheet', 'Summary_Centerline');
end

function saveVesselExcelData(vesselName, vesselNumber, pointOfInterest, data_struct, output_path)
    % Extract relevant data
    branchList = data_struct.branchList;
    Logical_branch = branchList(:, 4) == vesselNumber;
    indicesInBranch = find(Logical_branch);
    selectedPointIndex = indicesInBranch(pointOfInterest);

    % Define a 5-point range centered on the selected point
    index_range = max(selectedPointIndex - 2, 1):min(selectedPointIndex + 2, length(branchList));
    index_range(~Logical_branch(index_range)) = [];

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

    % Labels for time-averaged data
    Labels = (pointOfInterest - 2):(pointOfInterest + 2);
    Labels(Labels < 1) = 0;
    Labels(Labels > length(branchList)) = 0;
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
    selectedPointIndex = indicesInBranch(pointOfInterest);

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
    fig = figure('Visible', 'off');
    montage(FinalImage, 'Size', [subL, 4]);
    
    savename = sprintf('%s_Vessel_%d_Point_%d_Slicesview.jpg', vesselName, vesselNumber, pointOfInterest);
    saveas(fig, fullfile(output_path, savename));
    close(fig);
end
