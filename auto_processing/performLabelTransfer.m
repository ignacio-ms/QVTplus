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
    registeredImagePath = performSPMRegistration(QVT_mag, tofOrigResampled, tofOrigEICABCW);

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

    spm_write_vol(V, updatedBinarySegMatrix);
end

function [correspondenceDict, multiQVT] = generateCorrespondenceDict(folderPath, data_struct)
    % Generate a correspondence dictionary between eICAB labels and QVT labels
    %
    % Inputs:
    % - folderPath: Path to the folder containing QVT and multilabel files
    % - data_struct: Structure containing vessel and branch information
    %
    % Outputs:
    % - correspondenceDict: Dictionary mapping QVT labels to eICAB labels

    % Load multi-label segmentation
    multiQVT = spm_read_vols(spm_vol(fullfile(folderPath, 'multilabel_QVTseg.nii')));

    % Initialize correspondence dictionary
    correspondenceDict = struct();
    positions = data_struct.branchList(:, 1:3);  
    labels = data_struct.branchList(:, 4);

    % Find all vessel segments from QVT that belong to each of the eICAB labels
    for i = 1:size(positions, 1)
        x = round(positions(i, 1));
        y = round(positions(i, 2));
        z = round(positions(i, 3));

        if x > 0 && y > 0 && z > 0 && ...
           x <= size(multiQVT, 1) && y <= size(multiQVT, 2) && z <= size(multiQVT, 3)
            good_lab = multiQVT(x, y, z);
        else
            % Skip if the position is out of bounds in multiQVT
            continue;
        end

        if good_lab == 0
            continue; % Skip if no label found
        end

        fieldName = sprintf('good_lab_%d', good_lab);
        if ~isfield(correspondenceDict, fieldName)
            correspondenceDict.(fieldName) = [];
        end
        correspondenceDict.(fieldName) = [correspondenceDict.(fieldName); labels(i)];
    end

    % Resolve multiple QVT labels mapped to the same eICAB label
    labelOccurrences = correspondence_funcs('buildLabelOccurrences', correspondenceDict);
    correspondenceDict = correspondence_funcs('resolveLabelMappings', correspondenceDict, labelOccurrences);

    % Remove duplicate entries within each key
    correspondenceDict = correspondence_funcs('removeDuplicateEntries', correspondenceDict);

    % Rename keys to their actual vessel names
    segmentMapping = {
        'good_lab_1', 'LICA';
        'good_lab_2', 'RICA';
        'good_lab_7', 'LMCA';
        'good_lab_8', 'RMCA';
        'good_lab_5', 'LACA';
        'good_lab_6', 'RACA';
        'good_lab_3', 'BASI';
        'good_lab_4', 'COMM';
        'good_lab_9', 'COMM';
        'good_lab_10', 'COMM';
        'good_lab_100', 'SSSV';
        'good_lab_100', 'LTSV';
        'good_lab_100', 'RTSV';
        'good_lab_100', 'STRV'
    };

    for i = 1:size(segmentMapping, 1)
        goodLab = segmentMapping{i, 1};
        targetKey = segmentMapping{i, 2};
        if isfield(correspondenceDict, goodLab)
            if strcmp(targetKey, 'COMM')
                if ~isfield(correspondenceDict, 'COMM')
                    correspondenceDict.COMM = [];
                end
                correspondenceDict.COMM = [correspondenceDict.COMM; correspondenceDict.(goodLab)];
            else
                correspondenceDict.(targetKey) = correspondenceDict.(goodLab);
            end
        end
    end

    % Process PCA segments
    segmentLabels = {'good_lab_13', 'good_lab_14', 'good_lab_15', 'good_lab_16'};
    segments = cell(size(segmentLabels));
    maxLengths = zeros(size(segmentLabels));

    for i = 1:numel(segmentLabels)
        if isfield(correspondenceDict, segmentLabels{i})
            segments{i} = correspondenceDict.(segmentLabels{i});
            if isempty(segments{i})
                maxLengths(i) = 0;
            else
                maxLengths(i) = max(arrayfun(@(j) sum(data_struct.branchList(:, 4) == j), segments{i}));
            end
        else
            segments{i} = [];
            maxLengths(i) = 0;
        end
    end

    if maxLengths(3) > maxLengths(1)
        correspondenceDict.LPCA = [segments{1}; segments{3}];
    else
        correspondenceDict.LPCA = segments{1};
        correspondenceDict.LPC2 = segments{3};
    end

    if maxLengths(4) > maxLengths(2)
        correspondenceDict.RPCA = [segments{2}; segments{4}];
    else
        correspondenceDict.RPCA = segments{2};
        correspondenceDict.RPC2 = segments{4};
    end

    if isfield(correspondenceDict, 'COMM')
        correspondenceDict.COMM = unique(correspondenceDict.COMM);
    else
        correspondenceDict.COMM = [];
    end

    % Remove unused good_lab keys
    allFields = fieldnames(correspondenceDict);
    fieldsToRemove = allFields(startsWith(allFields, 'good_lab_'));
    correspondenceDict = rmfield(correspondenceDict, fieldsToRemove);

    % Remove empty fields
    fieldNames = fieldnames(correspondenceDict);
    for i = 1:numel(fieldNames)
        if isempty(correspondenceDict.(fieldNames{i}))
            correspondenceDict = rmfield(correspondenceDict, fieldNames{i});
        end
    end
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