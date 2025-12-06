function timeMIPcrossectionTR = computeTimeResolvedCD(data_struct, imageData, path_to_data, ~)
% computeTimeResolvedCD - Compute time-resolved Complex Difference cross-sections
%
% This function computes Complex Difference (CD) images for each cardiac frame
% and interpolates them to the cross-sectional planes defined by the vessel
% centerlines. This enables time-resolved visualization of CD in the GUI.
%
% Inputs:
%   data_struct  - Processed data structure containing branchList, matrix, etc.
%   imageData    - Image data structure containing MAG and potentially velocity data
%   path_to_data - Path to original data directory (for reloading velocity if needed)
%   output_path  - Output directory path
%
% Output:
%   timeMIPcrossectionTR - Time-resolved CD cross-sections [nPoints x width^2 x nframes]
%
% Usage:
%   Called from paramMap_auto.m after data loading and processing
%
% Notes:
%   - Uses time-averaged MAG for all frames (magnitude changes minimally over cardiac cycle)
%   - Computes CD for each frame using time-resolved velocity
%   - Interpolates each frame's CD to the cross-sectional planes
%   - If velocity data is not in imageData, it will be reloaded from NIfTI files

    fprintf('Computing time-resolved Complex Difference cross-sections...\n');
    
    % Extract parameters from data_struct
    branchList = data_struct.branchList;
    matrix = data_struct.matrix;
    nframes = data_struct.nframes;
    VENC = data_struct.VENC;
    r = data_struct.r;
    
    % Get MAG from imageData
    if isfield(imageData, 'MAG')
        MAG = imageData.MAG;
    else
        error('MAG not found in imageData');
    end
    
    % Recreate plane interpolation parameters (same as in paramMap_params_threshS)
    InterpVals = 4;
    width = r * InterpVals * 2 + 1;
    segments = length(branchList);
    
    [x_full, y_full, z_full, x, y, z, ~, ~] = create_planes(branchList, r, single(matrix), InterpVals, width);
    
    % Try to get velocity data from imageData, otherwise reload
    if isfield(imageData, 'v') && ~isempty(imageData.v)
        v = imageData.v;
        fprintf('  Using velocity data from imageData\n');
    else
        fprintf('  Velocity data not in imageData, reloading from NIfTI files...\n');
        % Reload velocity data (similar to loadNII_auto)
        v = reloadVelocityData(path_to_data, matrix, nframes);
        if isempty(v)
            error('Unable to reload velocity data. Please ensure NIfTI files are available.');
        end
    end
    
    % Initialize output array
    timeMIPcrossectionTR = zeros(segments, width^2, nframes, 'single');
    
    % Process each frame
    for frame = 1:nframes
        if mod(frame, 5) == 0 || frame == 1
            fprintf('  Processing frame %d/%d\n', frame, nframes);
        end
        
        % Extract velocity for current frame
        % v is [a x c x b x 3 x nframes], extract frame
        vFrame = squeeze(v(:,:,:,:,frame)); % [a x c x b x 3]
        
        % Compute time-resolved Complex Difference for this frame
        % Using calc_angio: angio = MAG .* sin((pi/2 * Vmag) / VENC)
        % where Vmag = sqrt(sum(vFrame.^2, 4))
        Vmag = sqrt(sum(vFrame.^2, 4)); % Speed image [a x c x b]
        
        % Cap velocity at VENC
        idx = find(Vmag > VENC);
        if ~isempty(idx)
            Vmag(idx) = VENC;
        end
        
        % Create complex-difference angiogram for this frame
        timeMIP_frame = MAG .* sin((pi/2 * Vmag) / VENC);
        
        % Interpolate to cross-sectional planes
        [timeMIPcrossectionTR(:,:,frame)] = interp_vol_to_planes(...
            timeMIP_frame, x, y, z, x_full, y_full, z_full, width, segments);
    end
    
    fprintf('  Completed time-resolved CD computation.\n');
end

function v = reloadVelocityData(path_to_data, ~, ~)
% reloadVelocityData - Reload time-resolved velocity data from NIfTI files
% This is a helper function to reload velocity if not in imageData
    
    v = [];
    
    % Check if 'scans' subfolder exists
    if isfolder(fullfile(path_to_data, 'scans'))
        base_dir = fullfile(path_to_data, 'scans');
    else
        base_dir = path_to_data;
    end
    
    % Find direction folders
    if isfolder(fullfile(base_dir, 'scans'))
        folders = dir(base_dir);
        folders(ismember({folders.name}, {'.', '..'})) = [];
        folders = {folders([folders.isdir]).name};
    else
        all_dirs = dir(base_dir);
        all_dirs = all_dirs([all_dirs.isdir]);
        folders = {all_dirs.name};
    end
    
    folderAP = folders(~cellfun('isempty', regexp(folders, 'AP')));
    if isempty(folderAP), return; end
    folderAP = folderAP{1};
    
    folderRL = folders(~cellfun('isempty', regexp(folders, 'RL')));
    if isempty(folderRL), return; end
    folderRL = folderRL{1};
    
    folderFH = folders(~cellfun('isempty', regexp(folders, 'FH')));
    if isempty(folderFH), return; end
    folderFH = folderFH{1};
    
    % Get paths
    ap_path = fullfile(base_dir, folderAP);
    rl_path = fullfile(base_dir, folderRL);
    fh_path = fullfile(base_dir, folderFH);
    
    if isfolder(fullfile(ap_path, 'NIFTI'))
        ap_path = fullfile(ap_path, 'NIFTI');
        rl_path = fullfile(rl_path, 'NIFTI');
        fh_path = fullfile(fh_path, 'NIFTI');
    end
    
    % Load phase volumes
    try
        vxvol = dir(fullfile(ap_path, '*_ph.nii.gz'));
        if isempty(vxvol)
            vxvol = dir(fullfile(ap_path, '*_ph.nii'));
        end
        if isempty(vxvol), return; end
        
        vxvol = spm_vol(fullfile(vxvol(1).folder, vxvol(1).name));
        vx = spm_read_vols(vxvol);
        
        vyvol = dir(fullfile(rl_path, '*_ph.nii.gz'));
        if isempty(vyvol)
            vyvol = dir(fullfile(rl_path, '*_ph.nii'));
        end
        if isempty(vyvol), return; end
        
        vyvol = spm_vol(fullfile(vyvol(1).folder, vyvol(1).name));
        vy = spm_read_vols(vyvol);
        
        vzvol = dir(fullfile(fh_path, '*_ph.nii.gz'));
        if isempty(vzvol)
            vzvol = dir(fullfile(fh_path, '*_ph.nii'));
        end
        if isempty(vzvol), return; end
        
        vzvol = spm_vol(fullfile(vzvol(1).folder, vzvol(1).name));
        vz = spm_read_vols(vzvol);
        
        % Reconstruct velocity array (same as loadNII_auto)
        [a,c,b,d] = size(vx);
        v = zeros([a,c,b,3,d],'single');
        v(:,:,:,2,:) = -squeeze(vx(:,:,:,:))*10;  % AP -> Y component (negated)
        v(:,:,:,1,:) = squeeze(vy(:,:,:,:))*10;   % RL -> X component
        v(:,:,:,3,:) = -squeeze(vz(:,:,:,:))*10; % FH -> Z component (negated)
        
    catch ME
        warning('computeTimeResolvedCD:ReloadFailed', ...
                'Failed to reload velocity data: %s', ME.message);
        v = [];
    end
end

