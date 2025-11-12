%% Decompress .nii.gz
folderPath = 'E:\PhD\data\processed\QVT\BMRI765891';

fileList = dir(fullfile(folderPath, '*_CW.nii.gz'));
compressedFilePath = fullfile(folderPath, fileList(1).name);
decompressedFilePath = strrep(compressedFilePath, '.gz', '');
gunzip(compressedFilePath);

V = spm_vol(decompressedFilePath);
data = spm_read_vols(V);
V.fname = decompressedFilePath; 
spm_write_vol(V, data);

%% save QVT data for SPM

folderPath = 'E:\PhD\data\processed\QVT\BMRI223186';
fileList = dir(fullfile(folderPath, 'qvtData_*.mat'));
load(fullfile(folderPath, fileList(1).name));

dataMatrix = imageData.Segmented;  

V = struct();
V.fname = 'QVT_seg.nii';      
V.dim = size(dataMatrix);             
V.dt = [spm_type('float32'), 0];      
V.mat = eye(4);                       
V.mat(1,1) = data_struct.VoxDims(1);                       
V.mat(2,2) = data_struct.VoxDims(2);                       
V.mat(3,3) = data_struct.VoxDims(3);                       

V = spm_write_vol(V, dataMatrix);

disp('NIfTI file saved successfully.');

%% preprocess eICAB (180 rotation x axis)
fileList = dir(fullfile(folderPath, '*_CW.nii'));
transformedImagePath = fullfile(folderPath, fileList(1).name);

V = spm_vol(transformedImagePath);
transformedData = spm_read_vols(V);

rotated = imrotate3(transformedData, 180, [1 0 0], 'nearest', 'crop');

[folder, baseName, ext] = fileparts(transformedImagePath);
newFileName = fullfile(folder, [baseName, '_rot', ext]);
V.fname = newFileName;      

V = spm_write_vol(V, rotated);

%% recenter images to get initialization for registration

img_path = 'C:\Users\u149879\Documents\GitHub\QVT_seg.nii';
V = spm_vol(img_path);
new_origin = (V.dim(1:3) + 1) / 2;
V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 
spm_get_space(img_path, V.mat);

img_path = newFileName;
V = spm_vol(img_path);
new_origin = (V.dim(1:3) + 1) / 2;
V.mat(1:3, 4) = V.mat(1:3, 1:3) * -new_origin'; 
spm_get_space(img_path, V.mat);

%% spm

spm('defaults', 'fmri');
spm_jobman('initcfg');

refImagePath = 'C:\Users\u149879\Documents\GitHub\QVT_seg.nii'; 
sourceImagePath = newFileName; 

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {refImagePath};  
matlabbatch{1}.spm.spatial.coreg.estimate.source = {sourceImagePath};  
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';  
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [10 4 2 1];       
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.01 0.01 0.01 0.001 0.001 0.001 0.01 0.01 0.01]; 
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [3 3];

spm_jobman('run', matlabbatch);
clear matlabbatch;

matlabbatch{1}.spm.spatial.coreg.write.ref = {refImagePath};  
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceImagePath};  
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;    
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];    
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;          
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r_';     

spm_jobman('run', matlabbatch);

%% load transformed eICAB seg
[folder, baseName, ext] = fileparts(newFileName);
newBaseName = ['r_', baseName];
newFileName = fullfile(folder, [newBaseName, ext]);

transformedImagePath = newFileName;  
V = spm_vol(transformedImagePath);
multiseg = spm_read_vols(V);

%% Get eICAB labels into QVT mask

[rows, cols, slices] = ind2sub(size(imageData.Segmented), find(imageData.Segmented == 1));
binaryPoints = [rows, cols, slices];
[multiRows, multiCols, multiSlices] = ind2sub(size(multiseg), find(multiseg > 0));
multiLabels = multiseg(multiseg > 0);  
multiLabelPoints = [multiRows, multiCols, multiSlices];

%% Expand the multilabels within a short radius

distanceThreshold = 5; % Threshold for initial k-NN search
k = 5; % Number of neighbors

[nearestIndices, distances] = knnsearch(multiLabelPoints, binaryPoints, 'K', k);

validNeighbors = arrayfun(@(i) label_transfer('processNeighbours',nearestIndices(i, :), distances(i, :), ...
                                               multiLabels, distanceThreshold, k), ...
                          1:size(nearestIndices, 1), 'UniformOutput', false);

updatedBinarySegMatrix = double(imageData.Segmented);
majorityLabels = label_transfer('assignMajorityLabels',validNeighbors, updatedBinarySegMatrix, rows, cols, slices);

linearIndices = sub2ind(size(updatedBinarySegMatrix), rows, cols, slices);
updatedBinarySegMatrix(linearIndices) = majorityLabels;

updatedBinarySegMatrix = label_transfer('expandLabels',updatedBinarySegMatrix, 100, 10);

%% save multilabel QVT mask to see in slicer

% V = struct();
% 
% V.fname = fullfile(folderPath, 'multilabel_QVTseg.nii');      
% V.dim = size(updatedBinarySegMatrix);             
% V.dt = [spm_type('float32'), 0];     
% V.mat = eye(4);                       
% V.mat(1,1) = data_struct.VoxDims(1);                      
% V.mat(2,2) = data_struct.VoxDims(2);                 
% V.mat(3,3) = data_struct.VoxDims(3);               
% 
% V = spm_write_vol(V, updatedBinarySegMatrix);

%%  find correspondence between eICAB labels and QVT labels
%   LICA: 1
%   RICA: 2
%   BA: 3
%   LACA: 5
%   RACA: 6
%   LMCA: 7
%   RMCA: 8
%   LPCA: 11 y 13 (15)
%   RPCA: 12 y 14 (16)
%   COMMS: 4, 9, 10
%   dont matter: 17, 18

folderPath = 'E:\PhD\data\processed\QVT\BMRI463775';
fileList = dir(fullfile(folderPath, 'qvtData_*.mat'));
load(fullfile(folderPath, fileList(1).name));

multiQVT = spm_read_vols(spm_vol(fullfile(folderPath,'multilabel_QVTseg.nii')));

correspondenceDict = struct();
positions = data_struct.branchList(:, 1:3);  
labels = data_struct.branchList(:, 4);

%find all the vessel segments from QVT that belong to each of the eICAB
%labels
for i = 1:size(positions, 1)

    x = round(positions(i, 1));
    y = round(positions(i, 2));
    z = round(positions(i, 3));
    
    if x > 0 && y > 0 && z > 0 && x <= size(multiQVT, 1) && y <= size(multiQVT, 2) && z <= size(multiQVT, 3)
        good_lab = multiQVT(x, y, z);
    else
        % Skip if the position is out of bounds in multiQVT
        continue;
    end
    
    if good_lab == 0 
        continue;
    end
    
    if ~isfield(correspondenceDict, sprintf('good_lab_%d', good_lab))
        correspondenceDict.(sprintf('good_lab_%d', good_lab)) = [];
    end
    
    correspondenceDict.(sprintf('good_lab_%d', good_lab)) = [correspondenceDict.(sprintf('good_lab_%d', good_lab)); labels(i)];
end

%Since QVT vessel labels might be mapped to several of the GT labels, keep
%the QVT labels only in the good_lab with the most occurrences.
labelOccurrences = correspondence_funcs('buildLabelOccurrences',correspondenceDict);

% Assign each label to the key with the most instances, except 'good_lab_100'
correspondenceDict = correspondence_funcs('resolveLabelMappings',correspondenceDict, labelOccurrences);

% Remove duplicate entries within each key
correspondenceDict = correspondence_funcs('removeDuplicateEntries',correspondenceDict);

%% Rename keys to their actual vessel name
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
    'good_lab_10', 'COMM'
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

segmentLabels = {'good_lab_13', 'good_lab_14', 'good_lab_15', 'good_lab_16'};
segments = cell(size(segmentLabels));
maxLengths = zeros(size(segmentLabels));

for i = 1:numel(segmentLabels)
    if isfield(correspondenceDict, segmentLabels{i})
        segments{i} = correspondenceDict.(segmentLabels{i});
        maxLengths(i) = max(arrayfun(@(j) sum(data_struct.branchList(:, 4) == j), segments{i}));
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

allFields = fieldnames(correspondenceDict);
fieldsToRemove = allFields(startsWith(allFields, 'good_lab_'));
correspondenceDict = rmfield(correspondenceDict, fieldsToRemove);

fieldNames = fieldnames(correspondenceDict);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if isempty(correspondenceDict.(fieldName))
        correspondenceDict = rmfield(correspondenceDict, fieldName);
    end
end

%% find ICA & BA cutoff slice (so the 3 vessels are measured in the same plane)

% Initialize variables
ica_slice = 0;
LICA = correspondenceDict.LICA;
RICA = correspondenceDict.RICA;
BA = correspondenceDict.BASI;

% Extract branch information for each vessel
info_LICA = find_LOCs('extractBranchInfo', data_struct, LICA);
info_RICA = find_LOCs('extractBranchInfo', data_struct, RICA);
info_BA = find_LOCs('extractBranchInfo', data_struct, BA);

% Round Z values
info_LICA(:, 3) = round(info_LICA(:, 3));
info_RICA(:, 3) = round(info_RICA(:, 3));
info_BA(:, 3) = round(info_BA(:, 3));

% Find unique Z values
unique_z_LICA = unique(info_LICA(:, 3));
unique_z_RICA = unique(info_RICA(:, 3));
unique_z_BA = unique(info_BA(:, 3));

% Find the largest Z value with one point for each vessel
max_single_z = find_LOCs('findMaxSingleZ', unique_z_LICA, unique_z_RICA, unique_z_BA, info_LICA, info_RICA, info_BA);

% Update ICA slice and extract locations
if max_single_z > -inf
    ica_slice = max_single_z;
    LICA_LOC = find_LOCs('extractLocation', info_LICA, max_single_z, 3);
    RICA_LOC = find_LOCs('extractLocation', info_RICA, max_single_z, 3);
    BA_LOC = find_LOCs('extractLocation', info_BA, max_single_z, 3);

    % Ensure minimum LOC Z offset of 3
    LICA_LOC = find_LOCs('ensureMinZOffset', info_LICA, LICA_LOC, 3);
    RICA_LOC = find_LOCs('ensureMinZOffset', info_RICA, RICA_LOC, 3);
    BA_LOC = find_LOCs('ensureMinZOffset', info_BA, BA_LOC, 3);
end


%% find the LOC value for each vessel. Ugly ugly code but it works

LOCs = struct();
vessel_labels = fieldnames(correspondenceDict);

for i = 1:numel(vessel_labels)
    keyName = vessel_labels{i};

    if ismember(keyName, {'LPCA', 'RPCA', 'LMCA', 'RMCA', 'RACA', 'LACA'})
        LOCs.(keyName) = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT);

    elseif strcmp(keyName, 'LICA')
        LOCs.(keyName) = [LICA_LOC(1, 4), LICA_LOC(1, 5)];
    elseif strcmp(keyName, 'RICA')
        LOCs.(keyName) = [RICA_LOC(1, 4), RICA_LOC(1, 5)];
    elseif strcmp(keyName, 'BASI')
        LOCs.(keyName) = [BA_LOC(1, 4), BA_LOC(1, 5)];

    elseif strcmp(keyName, 'COMM')
        LOCs.(keyName) = correspondenceDict.COMM;
    end
end

%% Keep either PCAs or PC2s

if isfield(correspondenceDict, 'RPC2') && isfield(LOCs, 'RPCA')
    if ismember(LOCs.('RPCA')(1), correspondenceDict.('RPC2'))
        correspondenceDict.('RPCA') = [correspondenceDict.('RPCA'); correspondenceDict.('RPC2')];
        correspondenceDict = rmfield(correspondenceDict, 'RPC2');
    else
        correspondenceDict = rmfield(correspondenceDict, 'RPC2');
    end
end

if isfield(correspondenceDict, 'LPC2') && isfield(LOCs, 'LPCA')
    if ismember(LOCs.('LPCA')(1), correspondenceDict.('LPC2'))
        correspondenceDict.('LPCA') = [correspondenceDict.('LPCA'); correspondenceDict.('LPC2')];
        correspondenceDict = rmfield(correspondenceDict, 'LPC2');
    else
        correspondenceDict = rmfield(correspondenceDict, 'LPC2');
    end
end

%% Save tagging into a .xlsx for QVT

vesselLabels = {
    'Left ICA', 'LICA';
    'Right ICA', 'RICA';
    'Basilar', 'BASI';
    'Left MCA', 'LMCA';
    'Right MCA', 'RMCA';
    'Left PCA', 'LPCA';
    'Right PCA', 'RPCA';
    'Left ACA', 'LACA';
    'Right ACA', 'RACA'
};

summaryData = cell(size(vesselLabels, 1), 7);
columnNames = {'Vessel Label', 'Centerline', 'Notes', 'Max Vel < 700 cm/s', ...
               'Mean Flow ml/s', 'Pulsatility Index', 'Branch Number'};

for i = 1:size(vesselLabels, 1)
    vesselLabel = vesselLabels{i, 1};
    locKey = vesselLabels{i, 2};

    summaryData{i, 1} = vesselLabel; % "Vessel Label"
    
    if isfield(LOCs, locKey)
        summaryData{i, 2} = LOCs.(locKey)(2); % "Centerline"
        summaryData{i, 7} = LOCs.(locKey)(1); % "Branch Number"
        rowIndex = find(data_struct.branchList(:, 4) == LOCs.(locKey)(1) & data_struct.branchList(:, 5) == LOCs.(locKey)(2));

        summaryData{i, 3} = ''; % "Notes"
        if data_struct.maxVel_val(rowIndex) < 600
            summaryData{i, 4} = 'YES'; % "Max Vel < 700 cm/s"
        else
            summaryData{i, 4} = 'NO';
        end
        summaryData{i, 5} = data_struct.flowPerHeartCycle_val(rowIndex); % "Mean Flow ml/s"
        summaryData{i, 6} = data_struct.PI_val(rowIndex); % "Pulsatility Index"
    else
        summaryData{i, 2} = NaN; 
        summaryData{i, 3} = NaN;
        summaryData{i, 4} = NaN;
        summaryData{i, 5} = NaN;
        summaryData{i, 6} = NaN;
        summaryData{i, 7} = NaN; 
    end
end

summaryTable = cell2table(summaryData, 'VariableNames', columnNames);
writetable(summaryTable, 'SummaryParamTools.xlsx');

%the actual export also creates more sheets in the xlsx with averaged and
%time resolved data for each vessel (the measuring point + 2 before + 2
%after, but i think this should be enough for now.
%Also, we're missing the Cavernous part of the ICAs (should be easy to
%include) as well as the VAs (not so easy to include) + venous system
%(would need to have a multilabel seg of the venous system).
%% Save tagging into the .csv file for QVT+

arteryKeyMap = struct( ...
    'L_ICA', 'LICA', ...
    'R_ICA', 'RICA', ...
    'L_MCA', 'LMCA', ...
    'R_MCA', 'RMCA', ...
    'L_ACA', 'LACA', ...
    'R_ACA', 'RACA', ...
    'L_PCA', 'LPCA', ...
    'R_PCA', 'RPCA', ...
    'BA', 'BASI', ...
    'Exclude', 'COMM' ...
);

arteryNames = {'L ICA', 'R ICA', 'L MCA', 'R MCA', 'L ACA', 'R ACA', 'L PCA', 'R PCA', 'BA', 'Exclude'};

arteryList = {};
labelList = {};
locList = {};

for i = 1:numel(arteryNames)
    artery = arteryNames{i};
    internalKey = arteryKeyMap.(strrep(artery, ' ', '_'));

    labelStr = ' ';
    locStr = ' ';

    if ~isempty(internalKey) && isfield(correspondenceDict, internalKey)
        labelData = correspondenceDict.(internalKey);
        labelStr = sprintf('[%s]', strjoin(string(labelData), ', '));
    end

    if ~isempty(internalKey) && isfield(LOCs, internalKey)
        locData = LOCs.(internalKey);
        if numel(locData) > 1
            locStr = sprintf('[%d, %d]', locData(1), locData(2));
        else
            locStr = sprintf('%d', locData(1));
        end
    end

    arteryList{end+1} = artery;
    labelList{end+1} = labelStr;
    locList{end+1} = locStr;
end

T = table(arteryList', labelList', locList', 'VariableNames', {'Artery', 'Label', 'Loc'});
writetable(T, 'artery_labels_auto.csv');