clear; clc;

% addpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus');
% addpath(genpath('/home/imarcoss/BioImaging/src/neuroimaging/qvt_plus'));
% addpath('/home/imarcoss/MATLAB Add-Ons/SWR/spm12-r7771');
% savepath;

% Define base paths
patient_id = 'PESA10758400';
% base_path = '/home/imarcoss/NetVolumes/Tierra/LAB_VF-ICH/LAB/MCC LAB/_IgnacioMarcos/LabVF/PESA-Brain/';
base_path = '/data_local/LabVF/PESA-Brain/DATA/Batch1_FP/';
use_eicab_whole_brain = true;

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
    % path_to_data = fullfile(base_path, 'DATA/Nifti', current_patient_id, '4DFlow')
    % eICAB_path = fullfile(base_path, 'RESULTS/eICAB', current_patient_id)
    path_to_data = fullfile(base_path, 'Nifti_Reorganized', current_patient_id, '4DFlow')
    eICAB_path = fullfile(base_path, 'eICAB', current_patient_id)
    % output_path = fullfile(base_path, 'RESULTS/QVTPlus', current_patient_id)
    % output_path = fullfile('/data_local/LabVF/PESA-Brain/RESULTS/QVTPlus_Refactored/', current_patient_id)
    output_path = fullfile('/data_local/LabVF/PESA-Brain/DATA/Batch1_FP/QVTPlus/', current_patient_id, '_FixedVel')
    mkdir(output_path);
    disp(['--------------------------------']);

    % Load data
    [data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

    % Perform label transfer and preprocessing
    [correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, path_to_data, data_struct, use_eicab_whole_brain);

    % Generate LOCs
    [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

    % Compute time-resolved Complex Difference cross-sections
    % This enables time-resolved CD visualization in the GUI when "Sync Images" is enabled
    try
        timeMIPcrossectionTR = computeTimeResolvedCD(data_struct, imageData, path_to_data, output_path);
        % Add to data_struct for saving
        data_struct.timeMIPcrossectionTR = timeMIPcrossectionTR;
        % Update the saved .mat file
        matInfo = dir(fullfile(output_path, 'qvtData_ISOfix_*.mat'));
        if ~isempty(matInfo)
            [~, newestIdx] = max([matInfo.datenum]);
            matFile = fullfile(output_path, matInfo(newestIdx).name);
            % Load existing data
            savedData = load(matFile);
            % Update data_struct
            savedData.data_struct.timeMIPcrossectionTR = timeMIPcrossectionTR;
            % Save back
            save(matFile, '-struct', 'savedData', '-v7.3');
            disp('Time-resolved CD data saved successfully.');
        end
    catch ME
        warning('computeTimeResolvedCD:Failed', ...
                'Failed to compute time-resolved CD: %s\nThis feature will not be available in the GUI.', ...
                ME.message);
    end

    % Save vessel-specific data automatically
    saveVesselData(LOCs, data_struct, output_path);

    %Save data for qvt+
    generateQVTplus(correspondenceDict, LOCs, output_path)
    disp(['Processing completed successfully for patient: ' current_patient_id]);
    
    % Save last output path for GUI auto-launch (if --gui flag is used)
    % Save to a temporary file in the output directory
    lastOutputFile = fullfile(output_path, '.last_output_path.txt');
    fid = fopen(lastOutputFile, 'w');
    if fid ~= -1
        fprintf(fid, '%s', output_path);
        fclose(fid);
    end
end
