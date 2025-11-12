function [correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, refImage, data_struct)
    % Perform label transfer between eICAB and QVT masks
    % Includes 180-degree flip, SPM registration, and label transfer
    % Inputs:
    % - eICAB_path: Path to eICAB data
    % - output_path: Directory for saving output files
    % - data_struct: QVT data structure containing vessel and segmentation information

    % Step 1: Decompress the .nii.gz file and recenter
    % decompressedFilePath = decompressAndRecenter(eICAB_path);
    
    % select eICAB orig_resampled file
    fileList = dir(fullfile(eICAB_path, '**/*_resampled.nii*'));
    if isempty(fileList)
        error('No eICAB orig_resampled file found in the specified path.');
    end

    tofOrigResampled = fullfile(fileList(1).folder, fileList(1).name);

    % if it is .gz, decompress it
    if contains(fileList(1).name, '.gz')
        compressedFilePath = fullfile(fileList(1).folder, fileList(1).name);
        gunzip(compressedFilePath, output_path);
        tofOrigResampled = fullfile(output_path, strrep(fileList(1).name, '.gz', ''));
    else
        % copy the file to the output path
        cp(fullfile(fileList(1).folder, fileList(1).name), output_path);
        tofOrigResampled = fullfile(output_path, fileList(1).name);
    end

    % select eICAB orig_eICAB_CW file
    fileList = dir(fullfile(eICAB_path, '**/*_eICAB_CW.nii*'));
    if isempty(fileList)
        error('No eICAB orig_eICAB_CW file found in the specified path.');
    end

    tofOrigEICABCW = fullfile(fileList(1).folder, fileList(1).name);

    % if it is .gz, decompress it
    if contains(fileList(1).name, '.gz')
        compressedFilePath = fullfile(fileList(1).folder, fileList(1).name);
        gunzip(compressedFilePath, output_path);
        tofOrigEICABCW = fullfile(output_path, strrep(fileList(1).name, '.gz', ''));
    else
        % copy the file to the output path
        copyfile(fullfile(fileList(1).folder, fileList(1).name), output_path);
        tofOrigEICABCW = fullfile(output_path, fileList(1).name);
    end

    
    [QVT_path, QVT_mag] = saveQVTseg(output_path, imageData, data_struct);

    %% Save QVT labels as volume

    branchList = data_struct.branchList;
    x = branchList(:,1);
    y = branchList(:,2);
    z = branchList(:,3);
    val = branchList(:,4);

    % Round coordinates to ensure they are voxel indices
    x = round(x);
    y = round(y);
    z = round(z);

    % Initialize volume using original segmented dimensions
    volumeSize = size(imageData.Segmented);
    volume = zeros(volumeSize, 'single');

    % Filter out-of-bounds coordinates (just in case)
    inBounds = x >= 1 & x <= volumeSize(1) & ...
               y >= 1 & y <= volumeSize(2) & ...
               z >= 1 & z <= volumeSize(3);
    if ~all(inBounds)
        warning('[Branch Mask Creation] Some coordinates are out of bounds and will be ignored.');
    end
    x = x(inBounds);
    y = y(inBounds);
    z = z(inBounds);
    val = val(inBounds);

    % Fill in volume
    % for i = 1:length(x)
    %     volume(x(i), y(i), z(i)) = val(i);
    % end
    % Linear indices for faster assignment
    linearIdx = sub2ind(volumeSize, x, y, z);
    volume(linearIdx) = val;

    % Create a SPM volume header
    V = struct();
    V.fname = [output_path '/branch_mask.nii'];
    V.dim = size(volume);                  % Image dimensions
    V.dt = [16, 0];                        % Data type: 16 = float32
    V.mat = eye(4);                        % Affine matrix: identity (voxel space)
    V.descrip = 'Branch mask from coordinates';

    % Set voxel dimensions
    V.mat(1,1) = data_struct.VoxDims(1);
    V.mat(2,2) = data_struct.VoxDims(2);
    V.mat(3,3) = data_struct.VoxDims(3);

    % Compute and set the new origin (To match QVT_seg and QVT_MAG saved origins)
    new_origin = (V.dim(1:3) + 1) / 2; % Center of the image
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 

    % Write the volume
    spm_write_vol(V, volume);

    %%
    

    % Step 2: Perform 180-degree flip along the X-axis
    % flippedFilePath = flipImage180(decompressedFilePath);

    % Get the AP magnitude image from the refImage folder
    fileList = dir(fullfile(refImage, '**/*AP*.nii*'));
    fileList = fileList(~contains({fileList.name}, '_ph'));
    if isempty(fileList)
        error('No AP magnitude image found in the specified path.');
    end

    % Generate a maximum intensity projection (MIP) image from the AP magnitude image
    AP_image = fullfile(fileList(1).folder, fileList(1).name);
    MIP_image = generateMIP(AP_image, output_path);

    % Step 3: Perform SPM registration
    registeredImagePath = performFSLRegistration(QVT_mag, tofOrigResampled, tofOrigEICABCW);

    % Step 4: Perform label transfer and create multi-label QVT segmentation
    updatedBinarySegMatrix = transferLabels(registeredImagePath, imageData);

    % Step 5: Save the resulting segmentation
    saveMultiLabelQVT(updatedBinarySegMatrix, output_path, data_struct);

    % Step 6: Generate correspondence dictionary
    [correspondenceDict, multiQVT] = generateCorrespondenceDict(output_path, data_struct);


    disp('Label transfer completed successfully.');
end

%% Subfunctions
function decompressedFilePath = decompressAndRecenter(folderPath)
    % Decompress and recenter the .nii.gz file
    fileList = dir(fullfile(folderPath, '*_CW.nii.gz'));
    if isempty(fileList)
        error('No compressed eICAB data found in the specified path.');
    end

    compressedFilePath = fullfile(folderPath, fileList(1).name);
    decompressedFilePath = strrep(compressedFilePath, '.gz', '');
    gunzip(compressedFilePath);

    V = spm_vol(decompressedFilePath);
    data = spm_read_vols(V);

    V = spm_vol(decompressedFilePath);
    new_origin = (V.dim(1:3) + 1) / 2;
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 
    spm_write_vol(V, data);
end

function [QVT_path, QVT_mag] = saveQVTseg(output_path, imageData, data_struct)

    dataMatrix = imageData.Segmented;
    V = struct();
    V.fname = fullfile(output_path, 'QVT_seg.nii');      
    V.dim = size(dataMatrix);             
    V.dt = [spm_type('float32'), 0];      
    V.mat = eye(4);                       

    % Set voxel dimensions
    V.mat(1,1) = data_struct.VoxDims(1);                       
    V.mat(2,2) = data_struct.VoxDims(2);                       
    V.mat(3,3) = data_struct.VoxDims(3);                         

    % Compute and set the new origin
    new_origin = (V.dim(1:3) + 1) / 2; % Center of the image
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 

    % Write the volume with the recentered transformation
    spm_write_vol(V, dataMatrix);

    % Output the path to the saved file
    QVT_path = V.fname;

    % Display confirmation
    disp(['NIfTI file saved and recentered successfully at: ', QVT_path]);

    % repeat for imageData.MAG
    dataMatrix = imageData.MAG;
    V = struct();
    V.fname = fullfile(output_path, 'QVT_MAG.nii');
    V.dim = size(dataMatrix);
    V.dt = [spm_type('float32'), 0];
    V.mat = eye(4);

    % Set voxel dimensions
    V.mat(1,1) = data_struct.VoxDims(1);
    V.mat(2,2) = data_struct.VoxDims(2);
    V.mat(3,3) = data_struct.VoxDims(3);

    % Compute and set the new origin
    new_origin = (V.dim(1:3) + 1) / 2; % Center of the image
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin';

    % Write the volume with the recentered transformation
    spm_write_vol(V, dataMatrix);

    % Output the path to the saved file
    QVT_mag = V.fname;

    % Display confirmation
    disp(['NIfTI file saved and recentered successfully at: ', QVT_mag]);

    % Additionally save Complex Difference image
    dataMatrix = imageData.CD;
    V = struct();
    V.fname = fullfile(output_path, 'QVT_CD.nii');
    V.dim = size(dataMatrix);
    V.dt = [spm_type('float32'), 0];
    V.mat = eye(4);

    % Set voxel dimensions
    V.mat(1,1) = data_struct.VoxDims(1);
    V.mat(2,2) = data_struct.VoxDims(2);
    V.mat(3,3) = data_struct.VoxDims(3);

    % Compute and set the new origin
    new_origin = (V.dim(1:3) + 1) / 2; % Center of the image
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin';

    % Write the volume with the recentered transformation
    spm_write_vol(V, dataMatrix);

    % Output the path to the saved file
    QVT_CD = V.fname;

    % Display confirmation
    disp(['NIfTI file saved and recentered successfully at: ', QVT_CD]);
end


function flippedFilePath = flipImage180(decompressedFilePath)
    % Perform 180-degree flip along the X-axis and recenter
    V = spm_vol(decompressedFilePath);
    data = spm_read_vols(V);

    % Apply 180-degree rotation along the X-axis
    flippedData = imrotate3(data, 180, [1 0 0], 'nearest', 'crop');

    % Update the transformation matrix for the flipped volume
    new_origin = (V.dim(1:3) + 1) / 2; % Compute new center
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; % Update the origin

    % Save the flipped image
    [folder, baseName, ext] = fileparts(decompressedFilePath);
    flippedFilePath = fullfile(folder, [baseName, ext]);
    V.fname = flippedFilePath;
    spm_write_vol(V, flippedData);
end


function registeredImagePath = performSPMRegistration(QVT_path, sourceImagePath, eICABImagePath)
    % first, flip sourceImage and eICABImage
    sourceImagePath = flipImage180(sourceImagePath);
    eICABImagePath = flipImage180(eICABImagePath);

    % Perform SPM-based registration to align eICAB and QVT masks
    spm_jobman('initcfg');

    % Estimate transformation
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {QVT_path};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sourceImagePath};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {eICABImagePath};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % % Apply transformation
    % matlabbatch{1}.spm.spatial.coreg.write.ref = {QVT_path};
    % matlabbatch{1}.spm.spatial.coreg.write.source = {eICAB_path};
    % matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    % matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    % matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    % matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r_';

    % spm_jobman('run', matlabbatch);

    % Return the registered image path
    [folder, baseName, ext] = fileparts(eICABImagePath);
    registeredImagePath = fullfile(folder, ['r_', baseName, ext]);
end

function registeredImagePath = performFSLRegistration(QVT_path, sourceImagePath, eICABImagePath)
    % first, flip sourceImage and eICABImage
    sourceImagePath = flipImage180(sourceImagePath);
    eICABImagePath = flipImage180(eICABImagePath);
    setenv("FSLOUTPUTTYPE", "NIFTI");

    % Perform FSL-based registration to align eICAB and QVT masks
    flirtCmd = sprintf('flirt -in %s -ref %s -omat %s -cost normmi -searchcost normmi -dof 6', ...
        sourceImagePath, QVT_path, fullfile(fileparts(eICABImagePath), 'transform.mat'));
    disp(flirtCmd)
    system(flirtCmd);

    % Define the registered image path
    [folder, baseName, ext] = fileparts(eICABImagePath);
    registeredImagePath = fullfile(folder, ['r_', baseName, ext]);
    
    % Apply the transformation to the eICAB image
    applyCmd = sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour', ...
        eICABImagePath, QVT_path, fullfile(fileparts(eICABImagePath), 'transform.mat'), ...
        registeredImagePath);
    system(applyCmd);
end

function updatedBinarySegMatrix = transferLabels(registeredImagePath, imageData)
    % Perform label transfer from eICAB to QVT
    V = spm_vol(registeredImagePath);
    multiseg = spm_read_vols(V);

    % Extract QVT binary points
    [rows, cols, slices] = ind2sub(size(imageData.Segmented), find(imageData.Segmented == 1));
    binaryPoints = [rows, cols, slices];

    % Extract multi-label points
    [multiRows, multiCols, multiSlices] = ind2sub(size(multiseg), find(multiseg > 0));
    multiLabels = multiseg(multiseg > 0);
    multiLabelPoints = [multiRows, multiCols, multiSlices];

    % Perform k-NN search for label transfer
    distanceThreshold = 5;
    k = 5;
    [nearestIndices, distances] = knnsearch(multiLabelPoints, binaryPoints, 'K', k);

    % Assign labels based on valid neighbors
    validNeighbors = arrayfun(@(i) label_transfer('processNeighbours', nearestIndices(i, :), distances(i, :), ...
                                                  multiLabels, distanceThreshold, k), ...
                              1:size(nearestIndices, 1), 'UniformOutput', false);

    % Assign majority labels
    updatedBinarySegMatrix = double(imageData.Segmented);
    majorityLabels = label_transfer('assignMajorityLabels', validNeighbors, updatedBinarySegMatrix, rows, cols, slices);

    % Apply the labels
    linearIndices = sub2ind(size(updatedBinarySegMatrix), rows, cols, slices);
    updatedBinarySegMatrix(linearIndices) = majorityLabels;

    % Expand the labels to nearby regions
    updatedBinarySegMatrix = label_transfer('expandLabels', updatedBinarySegMatrix, 100, 10);
end

function saveMultiLabelQVT(updatedBinarySegMatrix, output_path, data_struct)
    % Save the multi-label QVT segmentation as a .nii file
    V = struct();
    V.fname = fullfile(output_path, 'multilabel_QVTseg.nii');
    V.dim = size(updatedBinarySegMatrix);
    V.dt = [spm_type('float32'), 0];
    V.mat = eye(4);
    V.mat(1, 1) = data_struct.VoxDims(1);
    V.mat(2, 2) = data_struct.VoxDims(2);
    V.mat(3, 3) = data_struct.VoxDims(3);

    % Compute and set the new origin (To match QVT_seg and QVT_MAG saved origins)
    new_origin = (V.dim(1:3) + 1) / 2; % Center of the image
    V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 

    spm_write_vol(V, updatedBinarySegMatrix);
end

function [MIP_image] = generateMIP(AP_image, output_path)
    % Generate a maximum intensity projection (MIP) image from the AP magnitude image
    V = spm_vol(AP_image);
    data = spm_read_vols(V);

    % Compute the MIP image
    MIP_data = max(data, [], 4);

    % Save the MIP image
    [folder, baseName, ext] = fileparts(AP_image);

    if contains(baseName, '.nii')
        baseName = strrep(baseName, '.nii', '');
    end

    if isempty(ext)
        ext = '.nii';
    elseif strcmp(ext, '.gz')
        ext = '.nii';
    end

    MIP_image = fullfile(output_path, [baseName, '_MIP', ext]);

    W = V(1);
    W.fname = MIP_image;
    spm_write_vol(W, MIP_data);
end