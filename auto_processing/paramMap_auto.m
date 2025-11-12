clear; clc;

% addpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus');
% addpath(genpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus'));
% addpath('/home/imarcoss/MATLAB Add-Ons/SWR/spm12-r7771');
% savepath;

% Define paths
% path_to_data = 'E:\PhD\data\CNIC\BMRI198894\';
% output_path = 'C:\Users\u149879\Desktop\autoQVT';
patient_id = 'PESA10758400_A';
path_to_data = fullfile('/home/imarcoss/NetVolumes/Tierra/LAB_FSC/LAB/PERSONAL/imarcoss/PESA-Brain/DATA/Nifti', patient_id, '4DFlow');
% path_to_data = '/data_local/BioIT_Data/PESA-Brain/AlbertoExported/PESA10758400/scans';
% output_path = '/home/imarcoss/NetVolumes/Tierra/LAB_FSC/LAB/PERSONAL/imarcoss/PESA-Brain/Results/QVTPlus/PESA10758400/';
output_path = fullfile('/data_local/BioIT_Data/PESA-Brain/RESULTS', patient_id);
% output_path = '/data_local/BioIT_Data/PESA-Brain/RESULTS/PESA10758400_alberto';

% Load data
[data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

% Perform label transfer and preprocessing
% eICAB_path = 'E:\PhD\data\processed\QVT\BMRI198894';
% eICAB_path = fullfile('/home/imarcoss/NetVolumes/Tierra/LAB_FSC/LAB/PERSONAL/imarcoss/PESA-Brain/Results/Segmentation/eICAB', patient_id);
eICAB_path = fullfile('/home/imarcoss/NetVolumes/Tierra/LAB_FSC/LAB/PERSONAL/imarcoss/PESA-Brain/DATA/Nifti', patient_id, 'eICAB');
[correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, path_to_data, data_struct);

% Generate LOCs
[correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

% Save vessel-specific data automatically
saveVesselData(LOCs, data_struct, output_path);

%Save data for qvt+
generateQVTplus(correspondenceDict, LOCs, output_path)

disp('Processing completed successfully.');