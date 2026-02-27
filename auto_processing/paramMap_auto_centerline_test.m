function paramMap_auto_centerline_test(outputDir, varargin)
% paramMap_auto_centerline_test - Test centerlines from eICAB (arterial) + CD venous region
%
% Uses the same registration-skipping logic as paramMap_auto_erosion_test: loads existing
% .mat and does not re-run registration. Builds a new centerline by two paths:
%   - Arterial: eICAB mask in original TOF space (TOF_eICAB_WB.nii) -> centerline in TOF ->
%     transform points to 4D flow using transform.mat (FSL flirt). Float coords from
%     registration are converted to integer voxel indices (one point per 4D flow voxel per
%     segment) so tangents and cross-section planes are stable (TOF res >> 4D flow res).
%   - Venous: CD binary mask in venous search region (4D flow space) -> centerline in 4D flow.
% Merged arterial + venous centerline points form branchList_new -> full flow/area pipeline.
% If TOF_eICAB_WB.nii or transform.mat are missing, falls back to single combined mask
% (registered r_TOF_eICAB_CW + venous) in 4D flow.
% Results saved to outputDir/centerline_test/.
%
% Usage:
%   paramMap_auto_centerline_test(outputDir)
%   paramMap_auto_centerline_test(globalOutputPath)  % batch: process all PESA* subdirs
%
% Inputs:
%   outputDir - Single: path to one QVT+ output directory. Batch: path containing PESA* subdirs.
%
% Outputs (all under outputDir/centerline_test/):
%   - qvtData_ISOfix_centerline_*.mat
%   - segment_centerline_eICAB_venous.nii (combined arterial+venous mask)
%   - segment_centerline_eICAB_only.nii (arterial/eICAB mask only, for checking)
%   - segment_centerline_venous_only.nii (venous mask only, for checking)
%   - branch_mask.nii (centerline as volume, same as paramMap_auto)
%   - flow_PI_per_centerline_centerline_test.csv (and .xlsx): per-vessel centerline points with flow_ml_s, PI, area_val (cross-section), circularity at each point
%   - SummaryParamTool.xls, vessel slice views (*_Slicesview.jpg), LabelsQVT.csv (remade from new centerline)

    set(0, 'DefaultFigureVisible', 'off')

    if nargin < 1 || isempty(outputDir)
        outputDir = '/data_local/LabVF/PESA-Brain/RESULTS/Auxiliary/QVTPlus_ReLOCs/'; % PESA15792676/
    end
    outputDir = char(outputDir);
    if ~exist(outputDir, 'dir')
        error('paramMap_auto_centerline_test:DirectoryNotFound', 'Output directory not found: %s', outputDir);
    end  

    matInRoot = dir(fullfile(outputDir, 'qvtData_ISOfix_*.mat'));
    if ~isempty(matInRoot)
        patientDirs = {outputDir};
    else
        pesaDirs = dir(fullfile(outputDir, 'PESA*'));
        pesaDirs = pesaDirs([pesaDirs.isdir]);
        if isempty(pesaDirs)
            error('paramMap_auto_centerline_test:MissingData', ...
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
            warning('paramMap_auto_centerline_test:Skip', 'No qvtData_ISOfix_*.mat in %s. Skipping.', outputDir);
            continue;
        end

        centerlineTestDir = fullfile(outputDir, 'centerline_test');
        if ~exist(centerlineTestDir, 'dir')
            mkdir(centerlineTestDir);
        end

        [~, newestIdx] = max([matInfo.datenum]);
        matFile = fullfile(outputDir, matInfo(newestIdx).name);
        disp('Loading existing data (same registration-skipping as erosion test)...');
        loaded = load(matFile);
        data_struct = loaded.data_struct;
        imageData = loaded.imageData;
        if isfield(loaded, 'Vel_Time_Res')
            Vel_Time_Res = loaded.Vel_Time_Res;
        else
            Vel_Time_Res = [];
        end

        if ~isfield(imageData, 'Segmented') || ~isfield(imageData, 'V') || ~isfield(imageData, 'CD')
            error('paramMap_auto_centerline_test:MissingData', 'Need imageData.Segmented, .V, .CD in %s.', matFile);
        end

        segmentCD = logical(imageData.Segmented);
        vMean = imageData.V;
        timeMIP = imageData.CD;
        MAG = imageData.MAG;
        matrix = data_struct.matrix;
        nframes = data_struct.nframes;
        res = data_struct.res;
        if isfield(data_struct, 'VoxDims')
            sliceSpace = data_struct.VoxDims(3);
        else
            sliceSpace = res;
        end
        if isfield(imageData, 'v')
            v = imageData.v;
        else
            error('paramMap_auto_centerline_test:MissingData', 'Need imageData.v for time-resolved flow.');
        end

        %% Masks in 4D flow space (for saving NIfTIs and for venous centerline / fallback)
        sz = size(segmentCD);
        third_y = round(sz(2) / 3);
        venous_region = false(sz);
        venous_region(:, 1:third_y, :) = true;

        % Arterial mask in 4D flow (registered eICAB) - used for saving segment_centerline_eICAB_only.nii and fallback
        arterial_mask = false(sz);
        multilabelPath = fullfile(outputDir, 'r_TOF_eICAB_CW.nii');
        if exist(multilabelPath, 'file')
            multiQVT = spm_read_vols(spm_vol(multilabelPath));
            if ~isequal(size(multiQVT), sz)
                multiQVT = imresize3(multiQVT, sz);
            end
            arterial_mask = (multiQVT > 0) & ~venous_region;
        else
            arterial_mask = segmentCD & ~venous_region;
        end

        venous_mask = segmentCD & venous_region;
        segment_combined = arterial_mask | venous_mask;
        areaThresh = round(sum(segment_combined(:)) * 0.005);
        segment_combined = bwareaopen(segment_combined, max(1, areaThresh), 6);

        sortingCriteria = 3;
        spurLength = 8;
        tofPath = fullfile(outputDir, 'TOF_eICAB_WB.nii');
        xfmPath = fullfile(outputDir, 'transform.mat');
        useTwoPath = exist(tofPath, 'file') == 2 && exist(xfmPath, 'file') == 2;

        if useTwoPath
            %% Two-path: arterial centerline in TOF -> transform to 4D flow; venous centerline in 4D flow; merge
            % --- Arterial: TOF space ---
            V_tof = spm_vol(tofPath);
            mask_tof = spm_read_vols(V_tof) > 0;
            mask_tof = bwareaopen(mask_tof, max(1, round(sum(mask_tof(:)) * 0.005)), 6);
            vMean_tof = zeros([size(mask_tof), 3], 'single');  % no velocity in TOF; flow-direction sort arbitrary
            [~, ~, branchList_arterial_tof, ~] = feature_extraction(sortingCriteria, spurLength, vMean_tof, mask_tof);
            nArt = size(branchList_arterial_tof, 1);
            disp(['Centerline test: extracted ' num2str(nArt) ' arterial centerline points in TOF space.']);

            % Transform arterial points from TOF voxel to 4D flow voxel
            branchList_arterial_4df = branchList_arterial_tof;
            useFSL = false;
            refNiiPath = fullfile(outputDir, 'QVT_MAG.nii');
            if exist(refNiiPath, 'file')
                [err, ~] = system('which img2imgcoord');
                if err == 0
                    useFSL = true;
                end
            end
            if useFSL
                % FSL img2imgcoord: source coords (0-based voxel, one per line); ref = QVT_MAG (same as flirt)
                tmpIn = [tempname() '_src.txt'];
                tmpOut = [tempname() '_dest.txt'];
                fid = fopen(tmpIn, 'w');
                for p = 1:nArt
                    % branchList (dim1, dim2, dim3) = (y, x, z); FSL expects same dim order, 0-based
                    fprintf(fid, '%d %d %d\n', ...
                        round(branchList_arterial_tof(p, 1)) - 1, ...
                        round(branchList_arterial_tof(p, 2)) - 1, ...
                        round(branchList_arterial_tof(p, 3)) - 1);
                end
                fclose(fid);
                cmd = sprintf('img2imgcoord -src %s -dest %s -xfm %s -vox %s 2>/dev/null > %s', ...
                    tofPath, refNiiPath, xfmPath, tmpIn, tmpOut);
                [status, ~] = system(cmd);
                if status == 0 && exist(tmpOut, 'file')
                    try
                        % img2imgcoord may produce variable format; read numeric lines only
                        coords = readImg2imgcoordOutput(tmpOut);
                        if ~isempty(coords) && size(coords, 1) == nArt && size(coords, 2) >= 3
                            branchList_arterial_4df(:, 1:3) = coords(:, 1:3) + 1;  % 0-based -> 1-based
                        else
                            useFSL = false;
                        end
                    catch
                        useFSL = false;
                    end
                else
                    useFSL = false;
                end
                try delete(tmpIn); catch, end
                try delete(tmpOut); catch, end
            end
            if ~useFSL
                % Fallback: matrix in world (mm); use source image that was used in flirt (TOF_resampled) if available
                srcVolPath = tofPath;
                resampledPath = fullfile(outputDir, 'TOF_resampled.nii');
                if exist(resampledPath, 'file')
                    srcVolPath = resampledPath;
                end
                V_src = spm_vol(srcVolPath);
                if exist(resampledPath, 'file') && ~isequal(V_src.dim, V_tof.dim)
                    for p = 1:nArt
                        y1 = branchList_arterial_tof(p, 1); x1 = branchList_arterial_tof(p, 2); z1 = branchList_arterial_tof(p, 3);
                        w = V_tof.mat * [y1 - 1; x1 - 1; z1 - 1; 1];
                        v_src = V_src.mat \ w;
                        branchList_arterial_tof(p, 1:3) = [v_src(1)+1, v_src(2)+1, v_src(3)+1];
                    end
                end
                fid = fopen(xfmPath, 'r');
                if fid >= 0
                    xfm = fscanf(fid, '%f', [4 4])';
                    fclose(fid);
                    xfm = double(xfm);
                    if isfield(data_struct, 'OriginalAffine')
                        V_ref_mat = data_struct.OriginalAffine;
                    elseif isfield(imageData, 'OriginalAffine')
                        V_ref_mat = imageData.OriginalAffine;
                    else
                        V_ref_mat = eye(4);
                        V_ref_mat(1,1) = data_struct.VoxDims(1);
                        V_ref_mat(2,2) = data_struct.VoxDims(2);
                        V_ref_mat(3,3) = data_struct.VoxDims(3);
                        no = (matrix(1:3)' + 1) / 2;
                        V_ref_mat(1:3, 4) = V_ref_mat(1:3, 1:3) * (-no);
                    end
                    for p = 1:nArt
                        y1 = branchList_arterial_tof(p, 1); x1 = branchList_arterial_tof(p, 2); z1 = branchList_arterial_tof(p, 3);
                        src_vox = [y1 - 1; x1 - 1; z1 - 1; 1];
                        src_mm = V_src.mat * src_vox;
                        ref_mm = xfm * src_mm;
                        ref_vox = V_ref_mat \ [ref_mm(1:3); 1];
                        branchList_arterial_4df(p, 1:3) = [ref_vox(1) + 1, ref_vox(2) + 1, ref_vox(3) + 1];
                    end
                end
            end
            % Clip to 4D flow bounds (float coords from registration)
            inBounds = branchList_arterial_4df(:,1) >= 1 & branchList_arterial_4df(:,1) <= sz(1) & ...
                       branchList_arterial_4df(:,2) >= 1 & branchList_arterial_4df(:,2) <= sz(2) & ...
                       branchList_arterial_4df(:,3) >= 1 & branchList_arterial_4df(:,3) <= sz(3);
            branchList_arterial_4df = branchList_arterial_4df(inBounds, :);

            % Convert float registration coords to integer voxel (pixel) indices and remove
            % repeated points from downscaling. TOF is higher res (e.g. 0.25x0.25x0.5 mm) than
            % 4D flow (~0.77x0.77x0.8 mm), so many TOF points map to the same 4D flow voxel;
            % float coords and dense points cause unstable tangents and wrong cross-section
            % planes in create_planes / paramMap_GUI. We snap to voxel centres and keep one
            % point per voxel per segment so centerline is on the image grid.
            branchList_arterial_4df(:, 1:3) = round(branchList_arterial_4df(:, 1:3));
            branchList_arterial_4df(:, 1) = max(1, min(sz(1), branchList_arterial_4df(:, 1)));
            branchList_arterial_4df(:, 2) = max(1, min(sz(2), branchList_arterial_4df(:, 2)));
            branchList_arterial_4df(:, 3) = max(1, min(sz(3), branchList_arterial_4df(:, 3)));

            reduced = [];
            for segId = unique(branchList_arterial_4df(:, 4))'
                seg = branchList_arterial_4df(branchList_arterial_4df(:, 4) == segId, :);
                % Preserve order along segment (column 5 = original centerline index)
                seg = sortrows(seg, 5);
                vox = seg(:, 1:3);
                keep = true(size(seg, 1), 1);
                for k = 2:size(seg, 1)
                    if isequal(vox(k, :), vox(k-1, :))
                        keep(k) = false;
                    end
                end
                seg = seg(keep, :);
                seg(:, 5) = (1:size(seg, 1))';
                reduced = [reduced; seg]; %#ok<AGROW>
            end
            branchList_arterial_4df = reduced;
            nArtSeg = 0;
            if ~isempty(branchList_arterial_4df)
                nArtSeg = max(branchList_arterial_4df(:, 4));
            end
            disp(['Centerline test: arterial after float->int voxel + dedup: ' num2str(size(branchList_arterial_4df, 1)) ' points']);

            % --- Venous: 4D flow space ---
            venous_clean = bwareaopen(venous_mask, max(1, round(sum(venous_mask(:)) * 0.005)), 6);
            [~, ~, branchList_venous_4df, ~] = feature_extraction(sortingCriteria, spurLength, vMean, venous_clean);
            nVen = size(branchList_venous_4df, 1);
            disp(['Centerline test: extracted ' num2str(nVen) ' venous centerline points in 4D flow space.']);

            % Renumber venous segment_id so they follow arterial segments, then merge
            branchList_venous_4df(:, 4) = branchList_venous_4df(:, 4) + nArtSeg;
            branchList_new = [branchList_arterial_4df; branchList_venous_4df];

            % Drop segments with too few points so create_planes (d=2) never indexes out of bounds (needs >= 2*d+1 = 5 points per segment)
            MIN_POINTS_PER_SEGMENT = 5;
            segIds = branchList_new(:, 4);
            [uSeg, ~, ic] = unique(segIds);
            countPerSeg = accumarray(ic, 1);
            keepSeg = countPerSeg >= MIN_POINTS_PER_SEGMENT;
            keptSegIds = uSeg(keepSeg);
            rowKeep = ismember(segIds, keptSegIds);
            branchList_new = branchList_new(rowKeep, :);
            % Renumber segment IDs to be contiguous 1..N
            oldToNew = zeros(max(segIds), 1);
            oldToNew(keptSegIds) = 1:numel(keptSegIds);
            branchList_new(:, 4) = oldToNew(branchList_new(:, 4));
            nDropped = sum(~rowKeep);
            if nDropped > 0
                disp(['Centerline test: dropped ' num2str(nDropped) ' points in segments with < ' num2str(MIN_POINTS_PER_SEGMENT) ' points; ' num2str(size(branchList_new, 1)) ' points remaining.']);
            end

            disp(['Centerline test: merged ' num2str(size(branchList_new, 1)) ' points (arterial in TOF->4df + venous in 4df).']);
        else
            if ~exist(tofPath, 'file')
                warning('paramMap_auto_centerline_test:NoTOF', 'TOF_eICAB_WB.nii not found. Using combined mask in 4D flow.');
            end
            if ~exist(xfmPath, 'file')
                warning('paramMap_auto_centerline_test:NoXfm', 'transform.mat not found. Using combined mask in 4D flow.');
            end
            %% Fallback: single combined mask in 4D flow (registered eICAB + venous)
            [~, ~, branchList_new, ~] = feature_extraction(sortingCriteria, spurLength, vMean, segment_combined);
            disp(['Centerline test: extracted ' num2str(size(branchList_new, 1)) ' centerline points from combined mask (fallback).']);
        end

        % Smooth centerline per segment (same method as feature_extraction) for stable tangents
        smoothParameter = 0.3750;
        branchListSmooth = zeros(size(branchList_new));
        for n = 1:max(branchList_new(:, 4))
            branchActual = branchList_new(branchList_new(:, 4) == n, :);
            npts = size(branchActual, 1);
            xyz = [branchActual(:, 1)'; branchActual(:, 2)'; branchActual(:, 3)'];
            xyzp = zeros(3, npts);
            for k = 1:3
                pp = csaps(1:npts, xyz(k, :), smoothParameter);
                xyzp(k, :) = ppval(pp, 1:npts);
            end
            mask = branchList_new(:, 4) == n;
            branchListSmooth(mask, 1:3) = xyzp';
            branchListSmooth(mask, 4:5) = branchList_new(mask, 4:5);
        end
        branchList_new = branchListSmooth;

        %% Run full flow/area pipeline with new branchList (same logic as loadNII_auto -> paramMap_params_threshS)
        filetype = 'nii';
        back = zeros(size(vMean), 'single');
        BGPCdone = 1;
        directory = outputDir;
        Exseg = [];
        try
            [area_val, diam_val, flowPerHeartCycle_val, maxVel_val, PI_val, RI_val, flowPulsatile_val, ...
                velMean_val, VplanesAllx, VplanesAlly, VplanesAllz, r, timeMIPcrossection, segmentFull, ...
                vTimeFrameave, MAGcrossection, bnumMeanFlow, bnumStdvFlow, StdvFromMean, Planes, pixelSpace, segmentFullJS, maxVelFrame] ...
                = paramMap_params_threshS(filetype, branchList_new, matrix, timeMIP, vMean, back, ...
                BGPCdone, directory, nframes, res, MAG, v, sliceSpace, Exseg);
        catch ME
            disp(getReport(ME, 'extended'));
            warning('paramMap_auto_centerline_test:ParamMapFailed', 'paramMap_params_threshS failed: %s', ME.message);
            continue;
        end

        PIvel_val = [];
        if exist('maxVelFrame', 'var')
            mv_mean = mean(maxVelFrame, 2);
            mv_mean(mv_mean == 0) = eps;
            PIvel_val = (max(maxVelFrame, [], 2) - min(maxVelFrame, [], 2)) ./ mv_mean;
        end

        %% Build data_struct for centerline test
        data_struct_new = data_struct;
        data_struct_new.area_val = area_val;
        data_struct_new.diam_val = diam_val;
        data_struct_new.flowPerHeartCycle_val = flowPerHeartCycle_val;
        data_struct_new.maxVel_val = maxVel_val;
        data_struct_new.PI_val = PI_val;
        data_struct_new.RI_val = RI_val;
        data_struct_new.flowPulsatile_val = flowPulsatile_val;
        data_struct_new.velMean_val = velMean_val;
        data_struct_new.timeMIPcrossection = timeMIPcrossection;
        data_struct_new.segmentFull = segmentFull;
        data_struct_new.vTimeFrameave = vTimeFrameave;
        data_struct_new.MAGcrossection = MAGcrossection;
        data_struct_new.Planes = Planes;
        data_struct_new.bnumMeanFlow = bnumMeanFlow;
        data_struct_new.bnumStdvFlow = bnumStdvFlow;
        data_struct_new.StdvFromMean = StdvFromMean;
        data_struct_new.pixelSpace = pixelSpace;
        data_struct_new.segmentFullJS = segmentFullJS;
        data_struct_new.branchList = branchList_new;
        data_struct_new.directory = centerlineTestDir;
        if ~isempty(PIvel_val)
            data_struct_new.PIvel_val = PIvel_val;
        end

        % Update Vel_Time_Res with the new centerline's time-resolved velocity planes
        Vel_Time_Res.VplanesAllx = VplanesAllx;
        Vel_Time_Res.VplanesAlly = VplanesAlly;
        Vel_Time_Res.VplanesAllz = VplanesAllz;

        caseFilePath = fullfile(centerlineTestDir, ['qvtData_ISOfix_centerline_' datestr(now, 'ddmmmyyyy_HHMM') '.mat']);
        data_struct = data_struct_new;
        save(caseFilePath, 'data_struct', 'Vel_Time_Res', 'imageData', '-v7.3');
        disp(['Saved: ' caseFilePath]);

        %% Save segment masks and branch centerline volume
        V_seg = struct();
        V_seg.dim = sz;
        V_seg.dt = [spm_type('float32'), 0];
        if isfield(data_struct, 'OriginalAffine')
            V_seg.mat = data_struct.OriginalAffine;
        else
            V_seg.mat = eye(4);
            V_seg.mat(1,1) = data_struct.VoxDims(1);
            V_seg.mat(2,2) = data_struct.VoxDims(2);
            V_seg.mat(3,3) = data_struct.VoxDims(3);
            no = (sz(:) + 1) / 2;
            V_seg.mat(1:3, 4) = V_seg.mat(1:3, 1:3) * (-no);
        end
        V_seg.fname = fullfile(centerlineTestDir, 'segment_centerline_eICAB_venous.nii');
        spm_write_vol(V_seg, double(segment_combined));
        disp(['Saved combined mask: ' V_seg.fname]);
        V_seg.fname = fullfile(centerlineTestDir, 'segment_centerline_eICAB_only.nii');
        spm_write_vol(V_seg, double(arterial_mask));
        disp(['Saved arterial-only mask: ' V_seg.fname]);
        V_seg.fname = fullfile(centerlineTestDir, 'segment_centerline_venous_only.nii');
        spm_write_vol(V_seg, double(venous_mask));
        disp(['Saved venous-only mask: ' V_seg.fname]);

        % Save final centerline as mask (equivalent to branch_mask.nii in paramMap_auto)
        volumeSize = size(imageData.Segmented);
        branch_mask_vol = zeros(volumeSize, 'single');
        x = round(branchList_new(:, 1));
        y = round(branchList_new(:, 2));
        z = round(branchList_new(:, 3));
        val = single(branchList_new(:, 4));
        inBounds = x >= 1 & x <= volumeSize(1) & y >= 1 & y <= volumeSize(2) & z >= 1 & z <= volumeSize(3);
        x = x(inBounds); y = y(inBounds); z = z(inBounds); val = val(inBounds);
        if ~isempty(x)
            linearIdx = sub2ind(volumeSize, x, y, z);
            branch_mask_vol(linearIdx) = val;
        end
        V_bm = struct();
        V_bm.fname = fullfile(centerlineTestDir, 'branch_mask.nii');
        V_bm.dim = volumeSize;
        V_bm.dt = [spm_type('float32'), 0];
        if isfield(data_struct, 'OriginalAffine')
            V_bm.mat = data_struct.OriginalAffine;
        elseif isfield(imageData, 'OriginalAffine')
            V_bm.mat = imageData.OriginalAffine;
        else
            V_bm.mat = eye(4);
            V_bm.mat(1,1) = data_struct.VoxDims(1);
            V_bm.mat(2,2) = data_struct.VoxDims(2);
            V_bm.mat(3,3) = data_struct.VoxDims(3);
            no = (volumeSize(:) + 1) / 2;
            V_bm.mat(1:3, 4) = V_bm.mat(1:3, 1:3) * (-no);
        end
        spm_write_vol(V_bm, branch_mask_vol);
        disp(['Saved centerline mask: ' V_bm.fname]);

        % Copy multilabel so generateCorrespondenceDict can run on centerline_test
        multilabel_src = fullfile(outputDir, 'multilabel_QVTseg.nii');
        if exist(multilabel_src, 'file')
            copyfile(multilabel_src, fullfile(centerlineTestDir, 'multilabel_QVTseg.nii'));
        end

        % Remake and save all measurements outputs: SummaryParamTool.xls, slice views, LabelsQVT.csv, flow per centerline point
        try
            [correspondenceDict_cl, multiQVT_cl] = generateCorrespondenceDict(centerlineTestDir, data_struct_new);
            [correspondenceDict_cl, LOCs_cl] = generateLOCs(data_struct_new, correspondenceDict_cl, multiQVT_cl);
            % Normalize LOCs so second value is a valid centerline index (1..n) for that segment when
            % possible. Venous LOCs from resolveLongVenousSegment/resolveSSSVSTRV use (seg, row_index);
            % do not overwrite those with midIdx — convert row index to centerline index or keep as-is.
            branchList_new = data_struct_new.branchList;
            if ~isempty(LOCs_cl)
                locKeys = fieldnames(LOCs_cl);
                for kk = 1:numel(locKeys)
                    key = locKeys{kk};
                    v = LOCs_cl.(key);
                    if isempty(v) || ~isfinite(v(1))
                        continue;
                    end
                    seg = v(1);
                    segRows = find(branchList_new(:, 4) == seg);
                    nPt = numel(segRows);
                    if nPt < 1
                        continue;
                    end
                    midIdx = max(1, round(nPt / 2));
                    if numel(v) < 2
                        LOCs_cl.(key) = [seg, midIdx];
                        continue;
                    end
                    pt = v(2);
                    % (seg, centerline_pt_index): branchList has a row with segment_id=seg and col5=pt
                    rowMatch = find(branchList_new(:, 4) == seg & branchList_new(:, 5) == pt, 1);
                    if ~isempty(rowMatch)
                        continue;
                    end
                    % Venous: (seg, row_index) from resolveLongVenousSegment — pt is branchList row index
                    if pt >= 1 && pt <= size(branchList_new, 1) && branchList_new(pt, 4) == seg
                        clIdx = branchList_new(pt, 5);
                        LOCs_cl.(key) = [seg, clIdx];
                        continue;
                    end
                    % Truly invalid: use middle point
                    LOCs_cl.(key) = [seg, midIdx];
                end
            end
            if ~isempty(LOCs_cl)
                branchList_new = data_struct_new.branchList;
                locFields = fieldnames(LOCs_cl);
                vNamesCl = {};
                segIdsCl = [];
                clPointIdx = [];
                pointIdxFull = [];
                flowVal = [];
                piVal = [];
                areaVal = [];
                circularityVal = [];
                for f = 1:numel(locFields)
                    vName = locFields{f};
                    vesselInfo = LOCs_cl.(vName);
                    if isempty(vesselInfo) || numel(vesselInfo) < 2
                        continue;
                    end
                    vesselNumber = vesselInfo(1);
                    idxBranch = find(branchList_new(:, 4) == vesselNumber);
                    idxBranch = idxBranch(:)';
                    for k = 1:numel(idxBranch)
                        n = idxBranch(k);
                        vNamesCl{end+1} = vName;
                        segIdsCl(end+1) = vesselNumber;
                        clPointIdx(end+1) = k;
                        pointIdxFull(end+1) = n;
                        flowVal(end+1) = flowPerHeartCycle_val(n);
                        piVal(end+1) = PI_val(n);
                        areaVal(end+1) = area_val(n);
                        circularityVal(end+1) = diam_val(n);
                    end
                end
                if ~isempty(vNamesCl)
                    Tcl = table(vNamesCl', segIdsCl', clPointIdx', pointIdxFull', flowVal', piVal', areaVal', circularityVal', ...
                        'VariableNames', {'vessel_name', 'segment_id', 'centerline_point_index', 'point_index', 'flow_ml_s', 'PI', 'area_val', 'circularity'});
                    csvPath = fullfile(centerlineTestDir, 'flow_PI_per_centerline_centerline_test.csv');
                    writetable(Tcl, csvPath);
                    disp(['Saved flow and PI per centerline point: ' csvPath]);
                    xlsxPath = fullfile(centerlineTestDir, 'flow_PI_per_centerline_centerline_test.xlsx');
                    try
                        writetable(Tcl, xlsxPath);
                        disp(['Saved flow and PI per centerline point (Excel): ' xlsxPath]);
                    catch
                    end
                end
            end
            saveVesselData(LOCs_cl, data_struct_new, centerlineTestDir);
            generateQVTplus(correspondenceDict_cl, LOCs_cl, centerlineTestDir);
            disp('Saved vessel data: SummaryParamTool.xls, slice views, LabelsQVT.csv to centerline_test.');
        catch ME
            disp(getReport(ME, 'extended'));
            warning('paramMap_auto_centerline_test:VesselSaveFailed', ...
                'generateCorrespondenceDict/generateLOCs/saveVesselData failed: %s', ME.message);
        end

    end

    disp('Centerline test completed. Results in centerline_test/.');

end

function coords = readImg2imgcoordOutput(filePath)
    % Read img2imgcoord output: lines of space-separated numbers (skip non-numeric/empty).
    coords = [];
    if isempty(filePath) || ~exist(filePath, 'file')
        return;
    end
    fid = fopen(filePath, 'r');
    if fid < 0
        return;
    end
    rows = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        line = strtrim(line);
        if isempty(line)
            continue;
        end
        parts = regexp(line, '\s+', 'split');
        nums = [];
        for j = 1:numel(parts)
            v = str2double(parts{j});
            if ~isnan(v)
                nums(end+1) = v; %#ok<AGROW>
            end
        end
        if ~isempty(nums)
            rows{end+1} = nums; %#ok<AGROW>
        end
    end
    fclose(fid);
    if isempty(rows)
        return;
    end
    coords = vertcat(rows{:});
end
