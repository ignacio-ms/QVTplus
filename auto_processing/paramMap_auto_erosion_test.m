function paramMap_auto_erosion_test(outputDir, varargin)
% paramMap_auto_erosion_test - Same pipeline & LOCs, flow remeasured with eroded cross-section masks
%
% Reads all original outputs and remakes flow/area measurements at the SAME LOCs (same
% centerline points) as the original, but using an eroded version of the cross-section mask
% (segmentFull) at each point so that area (and thus flow) at that LOC is computed from a
% smaller lumen. Only the interpolated mask used for flow (segmentFull) is eroded; the
% original 3D binary mask is not modified. LOCs are unchanged. New measurements are saved
% in outputDir/erosion_test/ for comparison with the originals.
%
% Usage:
%   paramMap_auto_erosion_test(outputDir)
%   paramMap_auto_erosion_test(outputDir, 'erosionSize', 2)
%   paramMap_auto_erosion_test(globalOutputPath)  % batch: process all PESA* subdirs
%
% Inputs:
%   outputDir - Single patient: path to one QVT+ output directory (qvtData_ISOfix_*.mat, etc.).
%               Batch: path to a folder that contains one subfolder per patient (names starting
%               with PESA...); each subfolder must be a QVT+ output directory.
%
% Optional Name-Value Pairs:
%   'erosionSize' - Radius for erosion in pixels (default: 1). Cross-sections use disk strel;
%                   the 3D binary mask uses spherical strel. Integer >= 0.
%
% Outputs:
%   - outputDir/erosion_test/qvtData_ISOfix_erosion_*.mat  (new area, flow, PI, etc.; same LOCs)
%   - outputDir/erosion_test/flow_masks_used_r*.nii        (3D: union of cross-section masks used for flow)
%   - outputDir/erosion_test/flow_PI_per_centerline_r*.csv (per vessel: flow and PI at every centerline point)
%   - outputDir/erosion_test/weighted_image_segmentFull_r*.nii (3D: weighted image used for segmentFull)
%   - outputDir/erosion_test/LOC_pixel_count_erosion_r*.csv   (one row per vessel: vessel_name, segment_id, centerline_index, point_index, num_pixels_flow_mask)
%   - outputDir/erosion_test/per_vessel_LOC_crosssections_r*/  (PNG per vessel: cross-section at that vessel's LOC + red outline of segmentFull mask)
%   - outputDir/erosion_test/ Summary tables, vessel Excel, LabelsQVT.csv (from new measurements)
%
% Flow masks 3D volume: At each LOC, segmentFull is a 2D binary mask (81x81) in the tangent plane.
%   Plane points have 3D coords (x_full,y_full,z_full). Every plane point where the mask is 1 is
%   rasterized into the image grid (rounded to nearest voxel). The 3D volume is the union (value 1
%   where any LOC's flow mask projects).
%
% Weighted image 3D volume: Same formula as in segment_cross_section_thresh: weightIMAGE =
%   0.2*magSLICE + 0.8*cdSLICE + 0.2*velSLICE (each slice normalized to [0,1]). Computed per LOC
%   from interpolated CD, velocity, magnitude in the tangent plane; rasterized to 3D (max when
%   multiple LOCs hit the same voxel).

    %% Erosion size (variable at top; overridable by name-value pair)
    set(0, 'DefaultFigureVisible', 'off')
    EROSION_SIZE = 2;  % pixels (disk for cross-sections; sphere for 3D mask)

    if nargin < 1 || isempty(outputDir)
        % error('paramMap_auto_erosion_test:MissingOutputDir', 'outputDir is required.');
        % outputDir = '/data_local/LabVF/PESA-Brain/RESULTS/Auxiliary/QVTPlus_ReLOCs/';
        outputDir = '/data_local/LabVF/PESA-Brain/RESULTS/Auxiliary/QVTPlus_ReLOCs/PESA6943225/';
    end
    if nargin >= 2
        for i = 1:2:length(varargin)
            if strcmp(varargin{i}, 'erosionSize')
                EROSION_SIZE = max(0, round(double(varargin{i+1})));
            end
        end
    end

    outputDir = char(outputDir);
    if ~exist(outputDir, 'dir')
        error('paramMap_auto_erosion_test:DirectoryNotFound', 'Output directory not found: %s', outputDir);
    end

    % Single patient vs batch: if path contains qvtData_ISOfix_*.mat use as-is; else use PESA* subdirs
    matInRoot = dir(fullfile(outputDir, 'qvtData_ISOfix_*.mat'));
    if ~isempty(matInRoot)
        patientDirs = {outputDir};
    else
        pesaDirs = dir(fullfile(outputDir, 'PESA*'));
        pesaDirs = pesaDirs([pesaDirs.isdir]);
        if isempty(pesaDirs)
            error('paramMap_auto_erosion_test:MissingData', ...
                  'No qvtData_ISOfix_*.mat in %s and no PESA* subdirectories found.', outputDir);
        end
        patientDirs = arrayfun(@(x) fullfile(outputDir, x.name), pesaDirs, 'UniformOutput', false);
        disp(['Batch mode: found ' num2str(numel(patientDirs)) ' patient(s) (PESA*).']);
    end

    scriptDir = fileparts(mfilename('fullpath'));
    qvtRoot = fileparts(scriptDir);
    if isempty(strfind(path, qvtRoot))
        addpath(genpath(qvtRoot));
    end

    for iPatient = 1:numel(patientDirs)
        outputDir = patientDirs{iPatient};
        [~, patientName] = fileparts(outputDir);
        matInfo = dir(fullfile(outputDir, 'qvtData_ISOfix_*.mat'));
        if isempty(matInfo)
            warning('paramMap_auto_erosion_test:Skip', 'No qvtData_ISOfix_*.mat in %s. Skipping.', outputDir);
            continue;
        end

    erosionOutputDir = fullfile(outputDir, 'erosion_test');
    if ~exist(erosionOutputDir, 'dir')
        mkdir(erosionOutputDir);
    end

    disp('================================');
    disp(['Patient: ' patientName ' (' num2str(iPatient) '/' num2str(numel(patientDirs)) ')']);
    disp(['Output dir: ' outputDir]);
    disp(['Erosion test: ' erosionOutputDir]);
    disp(['Erosion size (pixels): ' num2str(EROSION_SIZE)]);
    disp('================================');

    [~, newestIdx] = max([matInfo.datenum]);
    matFile = fullfile(outputDir, matInfo(newestIdx).name);
    disp('Loading existing data...');
    loaded = load(matFile);
    data_struct = loaded.data_struct;
    Vel_Time_Res = loaded.Vel_Time_Res;
    imageData = loaded.imageData;
    LOCs = [];
    correspondenceDict = [];
    if isfield(loaded, 'LOCs')
        LOCs = loaded.LOCs;
    end
    if isfield(loaded, 'correspondenceDict')
        correspondenceDict = loaded.correspondenceDict;
    end
    % Use exact same LOCs as original: from .mat or from existing SummaryParamTool.xls (same structure: Centerline, Branch Number)
    if isempty(LOCs)
        summaryPath = fullfile(outputDir, 'SummaryParamTool.xls');
        if exist(summaryPath, 'file')
            % Load LOCs from original SummaryParamTool.xls (columns: Vessel Label, Centerline, Notes, Max Vel, Mean Flow, Pulsatility, Branch Number)
            try
                [~, ~, raw] = xlsread(summaryPath, 'Summary_Centerline');
                locKeys = {'LICA'; 'RICA'; 'BASI'; 'LMCA'; 'RMCA'; 'LPCA'; 'RPCA'; 'LACA'; 'RACA'; 'RCOMM'; 'LCOMM'; 'STRV'; 'SSSV'; 'LTSV'; 'RTSV'};
                LOCs = struct();
                for row = 2:min(size(raw, 1), 16)
                    branchNum = raw{row, 7};
                    centerline = raw{row, 2};
                    if row-1 <= numel(locKeys) && isnumeric(branchNum) && isnumeric(centerline) && ~isnan(branchNum) && ~isnan(centerline)
                        LOCs.(locKeys{row-1}) = [double(branchNum), double(centerline)];
                    end
                end
                disp('Using exact same LOCs as original (loaded from SummaryParamTool.xls: Branch Number and Centerline columns).');
            catch ME
                error('paramMap_auto_erosion_test:MissingLOCs', ...
                      'LOCs not in .mat and could not load from SummaryParamTool.xls: %s. Save LOCs in .mat or ensure SummaryParamTool.xls exists.', ME.message);
            end
        end
        if isempty(LOCs)
            error('paramMap_auto_erosion_test:MissingLOCs', ...
                  ['LOCs are not in the .mat file and SummaryParamTool.xls was not found. ' ...
                   'This script measures at the exact same LOCs as the original (Branch Number and Centerline). ' ...
                   'Run the main QVT+ pipeline to produce SummaryParamTool.xls, or use a .mat that has LOCs saved.']);
        end
    else
        disp('Using exact same LOCs as original measurements (loaded from .mat).');
    end

    if ~isfield(imageData, 'Segmented') || ~isfield(data_struct, 'segmentFull')
        error('paramMap_auto_erosion_test:MissingData', ...
              'Need imageData.Segmented and data_struct.segmentFull in %s.', matFile);
    end

    %% 1) Erode only the interpolated cross-section mask (segmentFull) used for flow measurement
    r = data_struct.r;
    InterpVals = 4;
    width = r * InterpVals * 2 + 1;
    segmentFull = data_struct.segmentFull;
    segments = size(segmentFull, 1);
    pixelSpace = data_struct.pixelSpace;
    res = data_struct.res;

    segmentFull_eroded = zeros(size(segmentFull), 'like', segmentFull);
    se_disk = strel('disk', EROSION_SIZE);
    for n = 1:segments
        seg = reshape(segmentFull(n, :), [width, width]);
        seg_eroded = imerode(logical(seg), se_disk);
        segmentFull_eroded(n, :) = seg_eroded(:);
    end

    %% 2) Recompute area_val and diam_val from eroded cross-sections (same formulas as segment_cross_section_thresh)
    area_val_eroded = zeros(segments, 1);
    diam_val_eroded = zeros(segments, 1);
    for n = 1:segments
        seg = reshape(segmentFull_eroded(n, :), [width, width]);
        dArea = (res/10) * (pixelSpace(n)/10);
        area_val_eroded(n, 1) = sum(seg(:)) * dArea * ((2*r+1)/(2*r*InterpVals+1))^2;
        if sum(seg(:)) > 0
            D = bwdist(~seg);
            Rin = max(D(:));
            [xLoc, yLoc] = find(bwperim(seg));
            if numel(xLoc) >= 2
                Dp = pdist2([xLoc, yLoc], [xLoc, yLoc]);
                Rout = max(Dp(:)) / 2;
                if Rout > 0
                    diam_val_eroded(n, 1) = Rin^2 / Rout^2;
                end
            end
        end
    end
    diam_val_eroded(isinf(diam_val_eroded) | isnan(diam_val_eroded)) = 0;

    %% 3) Recompute flow at same LOCs: same planes, velocity masked by eroded segmentFull
    branchList = data_struct.branchList;
    matrix = data_struct.matrix;
    nframes = data_struct.nframes; 
    v = imageData.v;
    vMean = imageData.V;

    [x_full, y_full, z_full, x, y, z, Tangent_V, ~] = create_planes(branchList, r, single(matrix), InterpVals, width);
    ROW = repmat((1:InterpVals:width)', [1 r*2+1]);
    COL = repmat(1:InterpVals*(width):(width)^2, [r*2+1 1]) - 1;
    idCOL = reshape(ROW + COL, [1 numel(ROW)]);

    flowPulsatile_val = zeros(segments, nframes);
    maxVelFrame = zeros(segments, nframes);
    velPulsatile_val = zeros(segments, nframes);  % used for velMean_val

    for j = 1:nframes
        vx = squeeze(v(:, :, :, 1, j));
        vy = squeeze(v(:, :, :, 2, j));
        vz = squeeze(v(:, :, :, 3, j));
        [v1] = interp_vol_to_planes(vx, x, y, z, x_full, y_full, z_full, width, segments);
        [v2] = interp_vol_to_planes(vy, x, y, z, x_full, y_full, z_full, width, segments);
        [v3] = interp_vol_to_planes(vz, x, y, z, x_full, y_full, z_full, width, segments);
        v1 = bsxfun(@times, v1, Tangent_V(:, 1));
        v2 = bsxfun(@times, v2, Tangent_V(:, 2));
        v3 = bsxfun(@times, v3, Tangent_V(:, 3));
        vTimeFrame = segmentFull_eroded .* (0.1 * (v1 + v2 + v3));
        denom = sum(vTimeFrame ~= 0, 2);
        vTimeFramerowMean = sum(vTimeFrame, 2) ./ (denom + (denom == 0));
        flowPulsatile_val(:, j) = abs(vTimeFramerowMean .* area_val_eroded);
        maxVelFrame(:, j) = max(vTimeFrame, [], 2);
        velPulsatile_val(:, j) = vTimeFramerowMean;
    end

    %% 4) Derived hemodynamic parameters (same formulas as paramMap_params_threshS)
    maxVel_val = max(maxVelFrame, [], 2);
    flowPerHeartCycle_val = sum(flowPulsatile_val, 2) / nframes;
    velMean_val = sum(velPulsatile_val, 2) / nframes;
    mf = mean(flowPulsatile_val, 2);
    mf(mf == 0) = eps;
    PI_val = abs(max(flowPulsatile_val, [], 2) - min(flowPulsatile_val, [], 2)) ./ mf;
    mx = max(flowPulsatile_val, [], 2);
    mx(mx == 0) = eps;
    RI_val = abs(max(flowPulsatile_val, [], 2) - min(flowPulsatile_val, [], 2)) ./ mx;

    bnumMeanFlow = zeros(max(branchList(:, 4)), 1);   % used in data_struct_new
    bnumStdvFlow = zeros(max(branchList(:, 4)), 1);   % used in data_struct_new
    for i = 1:max(branchList(:, 4))
        idx1 = branchList(:, 4) == i;
        bnumMeanFlow(i) = mean(flowPerHeartCycle_val(idx1));
        bnumStdvFlow(i) = std(flowPerHeartCycle_val(idx1));
    end

    StdvFromMean = flowPerHeartCycle_val;
    for n = 1:max(branchList(:, 4))
        IDbranch = find(branchList(:, 4) == n);
        len = length(IDbranch);
        L_ID = ones(1, len);
        if len >= 3
            L_ID(3:end) = L_ID(3:end) + (1:(len-2));
        end
        R_ID = 3:(len+2);
        R_ID((end-2):end) = len;
        for m = 1:len
            QV_meanflow_var = 1 - std(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m)))) / (abs(mean(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m))))) + eps);
            QV_area_var = 1 - (std(area_val_eroded(IDbranch(L_ID(m):R_ID(m)))) / (abs(mean(area_val_eroded(IDbranch(L_ID(m):R_ID(m))))) + eps));
            QV_circularity = mean(diam_val_eroded(IDbranch(L_ID(m):R_ID(m))));
            minmax_phase = zeros(1, nframes);
            for kk = 1:nframes
                flows_phase = flowPulsatile_val(IDbranch(L_ID(m):R_ID(m)), kk);
                minmax_phase(kk) = max(flows_phase) - min(flows_phase);
            end
            QV_tightness = 1 - mean(minmax_phase) / (abs(mean(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m))))) + eps);
            StdvFromMean(IDbranch(m)) = QV_meanflow_var + QV_area_var + QV_circularity + QV_tightness;
        end
    end

    mv_mean = mean(maxVelFrame, 2);
    mv_mean(mv_mean == 0) = eps;
    PIvel_val = (max(maxVelFrame, [], 2) - min(maxVelFrame, [], 2)) ./ mv_mean;

    %% 5) Build new data_struct (same LOCs, branchList, etc.; updated measurements only)
    data_struct_new = data_struct;
    data_struct_new.area_val = area_val_eroded;
    data_struct_new.diam_val = diam_val_eroded;
    data_struct_new.segmentFull = segmentFull_eroded;
    data_struct_new.flowPulsatile_val = flowPulsatile_val;
    data_struct_new.flowPerHeartCycle_val = flowPerHeartCycle_val;
    data_struct_new.maxVel_val = maxVel_val;
    data_struct_new.velMean_val = velMean_val;
    data_struct_new.PI_val = PI_val;
    data_struct_new.RI_val = RI_val;
    data_struct_new.bnumMeanFlow = bnumMeanFlow;
    data_struct_new.bnumStdvFlow = bnumStdvFlow;
    data_struct_new.StdvFromMean = StdvFromMean;
    data_struct_new.PIvel_val = PIvel_val;
    data_struct_new.directory = erosionOutputDir;

    %% 6) Save new .mat and copy multilabel so downstream find files in erosion_test
    caseFilePath = fullfile(erosionOutputDir, ['qvtData_ISOfix_erosion_' num2str(EROSION_SIZE) '_' datestr(now, 'ddmmmyyyy_HHMM') '.mat']);
    data_struct = data_struct_new;
    save(caseFilePath, 'data_struct', 'Vel_Time_Res', 'imageData', '-v7.3');
    disp(['Saved: ' caseFilePath]);

    multilabel_src = fullfile(outputDir, 'multilabel_QVTseg.nii');
    if exist(multilabel_src, 'file')
        copyfile(multilabel_src, fullfile(erosionOutputDir, 'multilabel_QVTseg.nii'));
    end

    %% 6b) Save 3D flow masks, 3D weighted image, CSV pixels per LOC, CSV flow and PI per centerline point per vessel
    volSize = data_struct.matrix(1:3);
    vol_flow_masks = zeros(volSize, 'single');
    vol_weighted = zeros(volSize, 'single');
    timeMIPcrossection = data_struct.timeMIPcrossection;
    MAGcrossection = data_struct.MAGcrossection;
    vTimeFrameave = data_struct.vTimeFrameave;
    weightIMS = [0.2 0.8 0.2];
    for n = 1:segments
        cdSLICE = reshape(timeMIPcrossection(n,:), [width, width]);
        cdSLICE = (cdSLICE - min(cdSLICE(:))) ./ (max(cdSLICE(:)) - min(cdSLICE(:)) + eps);
        velSLICE = reshape(vTimeFrameave(n,:), [width, width]);
        velSLICE = (velSLICE - min(velSLICE(:))) ./ (max(velSLICE(:)) - min(velSLICE(:)) + eps);
        magSLICE = reshape(MAGcrossection(n,:), [width, width]);
        magSLICE = (magSLICE - min(magSLICE(:))) ./ (max(magSLICE(:)) - min(magSLICE(:)) + eps);
        weightIMAGE = weightIMS(1)*magSLICE + weightIMS(2)*cdSLICE + weightIMS(3)*velSLICE;
        for k = 1:(width*width)
            ix = round(x_full(n, k)); iy = round(y_full(n, k)); iz = round(z_full(n, k));
            if ix >= 1 && ix <= volSize(1) && iy >= 1 && iy <= volSize(2) && iz >= 1 && iz <= volSize(3)
                if segmentFull_eroded(n, k) > 0
                    vol_flow_masks(ix, iy, iz) = 1;
                end
                vol_weighted(ix, iy, iz) = max(vol_weighted(ix, iy, iz), single(weightIMAGE(k)));
            end
        end
    end
    V_fm = struct();
    V_fm.fname = fullfile(erosionOutputDir, ['flow_masks_used_r' num2str(EROSION_SIZE) '.nii']);
    V_fm.dim = volSize; V_fm.dt = [spm_type('float32'), 0];
    if isfield(data_struct, 'OriginalAffine')
        V_fm.mat = data_struct.OriginalAffine;
    else
        V_fm.mat = eye(4);
        V_fm.mat(1,1) = data_struct.VoxDims(1); V_fm.mat(2,2) = data_struct.VoxDims(2); V_fm.mat(3,3) = data_struct.VoxDims(3);
        no = (volSize(:) + 1) / 2; V_fm.mat(1:3, 4) = V_fm.mat(1:3, 1:3) * (-no);
    end
    spm_write_vol(V_fm, vol_flow_masks);
    disp(['Saved flow masks 3D: ' V_fm.fname]);
    V_w = struct();
    V_w.fname = fullfile(erosionOutputDir, ['weighted_image_segmentFull_r' num2str(EROSION_SIZE) '.nii']);
    V_w.dim = volSize; V_w.dt = [spm_type('float32'), 0]; V_w.mat = V_fm.mat;
    spm_write_vol(V_w, vol_weighted);
    disp(['Saved weighted image 3D: ' V_w.fname]);
    num_pixels_per_loc = sum(segmentFull_eroded > 0, 2);

    % CSV: for each vessel, flow and PI at every centerline point (not only the LOC)
    vNamesCl = {};
    segIdsCl = [];
    clPointIdx = [];
    pointIdxFull = [];
    flowVal = [];
    piVal = [];
    locFields = fieldnames(LOCs);
    for f = 1:numel(locFields)
        vName = locFields{f};
        vesselInfo = LOCs.(vName);
        if isempty(vesselInfo) || numel(vesselInfo) < 2
            continue;
        end
        vesselNumber = vesselInfo(1);
        idxBranch = find(branchList(:, 4) == vesselNumber);
        idxBranch = idxBranch(:)';
        for k = 1:numel(idxBranch)
            n = idxBranch(k);
            vNamesCl{end+1} = vName;
            segIdsCl(end+1) = vesselNumber;
            clPointIdx(end+1) = k;
            pointIdxFull(end+1) = n;
            flowVal(end+1) = flowPerHeartCycle_val(n);
            piVal(end+1) = PI_val(n);
        end
    end
    if ~isempty(vNamesCl)
        Tcl = table(vNamesCl', segIdsCl', clPointIdx', pointIdxFull', flowVal', piVal', ...
            'VariableNames', {'vessel_name', 'segment_id', 'centerline_point_index', 'point_index', 'flow_ml_s', 'PI'});
        writetable(Tcl, fullfile(erosionOutputDir, ['flow_PI_per_centerline_r' num2str(EROSION_SIZE) '.csv']));
        disp(['Saved flow and PI per centerline point: flow_PI_per_centerline_r' num2str(EROSION_SIZE) '.csv']);
    end

    %% 7) Save vessel tables and LabelsQVT from new measurements (exact same LOCs as original)
    if isempty(correspondenceDict)
        try
            [correspondenceDict, ~] = generateCorrespondenceDict(erosionOutputDir, data_struct_new);
        catch
            warning('paramMap_auto_erosion_test:NoCorrespondence', ...
                    'Could not build correspondenceDict. Some vessel operations may be limited.');
        end
    end
    if ~isempty(LOCs)
        try
            saveVesselData(LOCs, data_struct_new, erosionOutputDir);
            if ~isempty(correspondenceDict)
                generateQVTplus(correspondenceDict, LOCs, erosionOutputDir);
            end
            disp('Vessel data and QVT+ labels (eroded measurements) saved to erosion_test.');
        catch ME
            warning('paramMap_auto_erosion_test:SaveFailed', 'saveVesselData/generateQVTplus failed: %s', ME.message);
        end
        % LOC pixel count CSV: one row per vessel (final LOC only)
        locFields = fieldnames(LOCs);
        vNames = {};
        segIds = [];
        clIdxs = [];
        ptIdxs = [];
        npix = [];
        for f = 1:numel(locFields)
            vName = locFields{f};
            vesselInfo = LOCs.(vName);
            if isempty(vesselInfo) || numel(vesselInfo) < 2
                continue;
            end
            vesselNumber = vesselInfo(1);
            pointOfInterest = vesselInfo(2);
            n = find(branchList(:, 4) == vesselNumber & branchList(:, 5) == pointOfInterest, 1);
            if isempty(n)
                continue;
            end
            vNames{end+1} = vName;
            segIds(end+1) = vesselNumber;
            clIdxs(end+1) = pointOfInterest;
            ptIdxs(end+1) = n;
            npix(end+1) = num_pixels_per_loc(n);
        end
        if ~isempty(vNames)
            T = table(vNames', segIds', clIdxs', ptIdxs', npix', ...
                'VariableNames', {'vessel_name', 'segment_id', 'centerline_index', 'point_index', 'num_pixels_flow_mask'});
            writetable(T, fullfile(erosionOutputDir, ['LOC_pixel_count_erosion_r' num2str(EROSION_SIZE) '.csv']));
            disp('Saved LOC pixel count CSV (one row per vessel).');
        end
        % Save cross-section image with segmentFull outline only for the final LOC per vessel
        perLocDir = fullfile(erosionOutputDir, ['per_vessel_LOC_crosssections_r' num2str(EROSION_SIZE)]);
        if ~exist(perLocDir, 'dir')
            mkdir(perLocDir);
        end
        timeMIPcrossection = data_struct.timeMIPcrossection;
        MAGcrossection = data_struct.MAGcrossection;
        vTimeFrameave = data_struct.vTimeFrameave;
        weightIMS = [0.2 0.8 0.2];
        locFields = fieldnames(LOCs);
        for f = 1:numel(locFields)
            vName = locFields{f};
            vesselInfo = LOCs.(vName);
            if isempty(vesselInfo) || numel(vesselInfo) < 2
                continue;
            end
            vesselNumber = vesselInfo(1);
            pointOfInterest = vesselInfo(2);
            % Row index in branchList for this LOC: segment_id == vesselNumber, centerline_index == pointOfInterest
            n = find(branchList(:, 4) == vesselNumber & branchList(:, 5) == pointOfInterest, 1);
            if isempty(n)
                continue;
            end
            cdSLICE = reshape(timeMIPcrossection(n,:), [width, width]);
            cdSLICE = (cdSLICE - min(cdSLICE(:))) ./ (max(cdSLICE(:)) - min(cdSLICE(:)) + eps);
            velSLICE = reshape(vTimeFrameave(n,:), [width, width]);
            velSLICE = (velSLICE - min(velSLICE(:))) ./ (max(velSLICE(:)) - min(velSLICE(:)) + eps);
            magSLICE = reshape(MAGcrossection(n,:), [width, width]);
            magSLICE = (magSLICE - min(magSLICE(:))) ./ (max(magSLICE(:)) - min(magSLICE(:)) + eps);
            weightIMAGE = weightIMS(1)*magSLICE + weightIMS(2)*cdSLICE + weightIMS(3)*velSLICE;
            bg = mat2gray(weightIMAGE);
            mask_2d = reshape(segmentFull_eroded(n,:), [width, width]);
            outline = bwperim(mask_2d, 8);
            R = bg; G = bg; B = bg;
            R(outline) = 1; G(outline) = 0; B(outline) = 0;
            imwrite(cat(3, R, G, B), fullfile(perLocDir, [vName '.png']));
        end
        disp(['Saved per-vessel LOC cross-section images (mask outline): ' perLocDir]);
    end

    end  % for iPatient

    disp('================================');
    disp('Erosion test completed. Compare erosion_test/ with original outputs.');
    disp('================================');
end
