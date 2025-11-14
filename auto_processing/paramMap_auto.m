clear; clc;

% addpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus');
% addpath(genpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus'));
% addpath('/home/imarcoss/MATLAB Add-Ons/SWR/spm12-r7771');
% savepath;

% Define base paths
patient_id = 'PESA10758400_A';
base_path = '/home/imarcoss/NetVolumes/Tierra/LAB_VF-ICH/LAB/MCC LAB/_IgnacioMarcos/LabVF/PESA-Brain/';

if strcmp(patient_id, 'all')
    dir_listing = dir(fullfile(base_path, 'DATA/Nifti'));
    dir_listing = dir_listing([dir_listing.isdir]);
    patient_ids = {dir_listing.name};
    patient_ids = patient_ids(~ismember(patient_ids, {'.', '..'}));
else
    patient_ids = {patient_id};
end

for idx = 1:numel(patient_ids)
    current_patient_id = patient_ids{idx};
    disp(['--------------------------------']);
    disp(['Processing patient: ' current_patient_id]);
    path_to_data = fullfile(base_path, 'DATA/Nifti', current_patient_id, '4DFlow')
    eICAB_path = fullfile(base_path, 'RESULTS/eICAB', current_patient_id)
    % output_path = fullfile(base_path, 'RESULTS/QVTPlus', current_patient_id)
    output_path = fullfile('/data_local/LabVF/PESA-Brain/RESULTS/', current_patient_id)
    mkdir(output_path);
    disp(['--------------------------------']);

    % Load data
    [data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

    % Perform label transfer and preprocessing
    [correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, path_to_data, data_struct);

    % Generate LOCs
    [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

    % Save vessel-specific data automatically
    saveVesselData(LOCs, data_struct, output_path);

    %Save data for qvt+
    generateQVTplus(correspondenceDict, LOCs, output_path)
    disp(['Processing completed successfully for patient: ' current_patient_id]);
end
