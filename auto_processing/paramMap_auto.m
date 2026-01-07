function paramMap_auto(varargin)
% paramMap_auto - Automated QVT+ processing pipeline
%
% Usage:
%   paramMap_auto()          % Uses default paths 
%   paramMap_auto('--help')  % Display help message
%
%   Single patient mode:
%   paramMap_auto(path_to_data, eICAB_path, output_path)
%   paramMap_auto(path_to_data, eICAB_path, output_path, use_eicab_whole_brain)
%
%   Batch mode:
%   paramMap_auto('--batch', base_path, patient_id, use_eicab_whole_brain, skip_existing)
%   paramMap_auto('--batch', base_path, 'all', true, true)
%
%   Batch mode folder structure:
%   base_path: Input
%   4DFlow: base_path/DATA/Nifti/patient_id/4DFlow
%   eICAB: base_path/RESULTS/eICAB/patient_id
%   output_path: base_path/RESULTS/QVTPlus/patient_id
%
% Arguments:
%   path_to_data - Path to 4DFlow data directory (single patient mode)
%   eICAB_path - Path to eICAB results directory (single patient mode)
%   output_path - Path to output directory (single patient mode)
%   base_path - Base path for batch processing (batch mode)
%   patient_id - Patient ID or 'all' for batch processing (batch mode)
%   use_eicab_whole_brain - Boolean flag for eICAB whole brain (default: true)
%   skip_existing - Boolean flag to skip existing outputs (default: true)
%
% Examples:
%   % Single patient
%   paramMap_auto('/data/sub-001/4DFlow', '/data/sub-001/eICAB', '/results/sub-001')
%
%   % Batch processing all patients
%   paramMap_auto('--batch', '/data_local/LabVF/PESA-Brain', 'all', true, true)
%
%   % Single patient with flags
%   paramMap_auto('/data/sub-001/4DFlow', '/data/sub-001/eICAB', '/results/sub-001', false)

% clear; clc;

% Parse command-line arguments
if nargin == 0
    % No arguments - use defaults 
    patient_id = 'all';
    base_path = '/data_local/LabVF/PESA-Brain';
    use_eicab_whole_brain = false;
    skip_existing = false;
    batch_mode = true;
    single_patient_mode = false;
elseif nargin == 1 && (strcmp(varargin{1}, '--help') || strcmp(varargin{1}, '-h') || strcmp(varargin{1}, '-?'))
    % Display help
    help('paramMap_auto');
    return;
elseif nargin >= 1 && strcmp(varargin{1}, '--batch')
    % Batch mode: --batch base_path patient_id use_eicab_whole_brain skip_existing
    if nargin < 3
        error('Batch mode requires at least: --batch base_path patient_id');
    end
    base_path = varargin{2};
    patient_id = varargin{3};
    if nargin >= 4
        use_eicab_whole_brain = logical(str2double(varargin{4}));
    else
        use_eicab_whole_brain = true;
    end
    if nargin >= 5
        skip_existing = logical(str2double(varargin{5}));
    else
        skip_existing = true;
    end
    batch_mode = true;
    single_patient_mode = false;
elseif nargin >= 3
    % Single patient mode: path_to_data eICAB_path output_path [use_eicab_whole_brain]
    path_to_data = varargin{1};
    eICAB_path = varargin{2};
    output_path = varargin{3};
    if nargin >= 4
        use_eicab_whole_brain = logical(str2double(varargin{4}));
    else
        use_eicab_whole_brain = false;
    end
    batch_mode = false;
    single_patient_mode = true;
    skip_existing = false; % Not applicable for single patient mode
else
    error('Invalid arguments. Use --help for usage information.');
end

% Process based on mode
if single_patient_mode
    % Single patient mode
    patient_ids = {'single_patient'};
    total_patients = 1;
    correct_patients = 0;
    failed_patients = 0;
    
    try
        current_patient_id = 'single_patient';
        disp(['--------------------------------']);
        disp(['Processing single patient']);
        disp(['Path to data: ' path_to_data]);
        disp(['eICAB path: ' eICAB_path]);
        disp(['Output path: ' output_path]);
        disp(['--------------------------------']);
        
        % Validate paths
        if ~exist(path_to_data, 'dir')
            error('Path to data does not exist: %s', path_to_data);
        end
        if ~exist(eICAB_path, 'dir')
            error('eICAB path does not exist: %s', eICAB_path);
        end
        if ~exist(output_path, 'dir')
            mkdir(output_path);
            disp(['Created output directory: ' output_path]);
        end
        
        % Process single patient
        processSinglePatient(path_to_data, eICAB_path, output_path, use_eicab_whole_brain);
        
        % Save last output path for GUI auto-launch
        lastOutputFile = fullfile(output_path, '.last_output_path.txt');
        fid = fopen(lastOutputFile, 'w');
        if fid ~= -1
            fprintf(fid, '%s', output_path);
            fclose(fid);
        end
        
        correct_patients = 1;
        disp(['Processing completed successfully']);
    catch ME
        disp(['Error processing patient']);
        disp(ME.message);
        if ~isempty(ME.stack)
            for k = 1:length(ME.stack)
                fprintf('  File: %s, Line: %d, Function: %s\n', ...
                    ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
            end
        end
        failed_patients = 1;
    end
else
    % Batch mode - process multiple patients
    if strcmp(patient_id, 'all')
        dir_listing = dir(fullfile(base_path, '/DATA/Nifti'));
        dir_listing = dir_listing([dir_listing.isdir]);
        patient_ids = {dir_listing.name};
        patient_ids = patient_ids(~ismember(patient_ids, {'.', '..'}));
    else
        patient_ids = {patient_id};
    end
    
    total_patients = numel(patient_ids);
    correct_patients = 0;
    failed_patients = 0;
    
    for idx = 1:numel(patient_ids)
        try
            current_patient_id = patient_ids{idx};
            disp(['--------------------------------']);
            disp(['Processing patient: ' current_patient_id]);
            path_to_data = fullfile(base_path, '/DATA/Nifti', current_patient_id, '4DFlow');
            eICAB_path = fullfile(base_path, 'RESULTS/eICAB', current_patient_id);
            output_path = fullfile(base_path, '/RESULTS/QVTPlus/', current_patient_id);
            if skip_existing && exist(output_path, 'dir')
                disp(['Skipping patient: ' current_patient_id ' because it already exists']);
                continue;
            end
            disp(['--------------------------------']);
            
            % Process single patient
            processSinglePatient(path_to_data, eICAB_path, output_path, use_eicab_whole_brain);
            
            % Save last output path for GUI auto-launch
            lastOutputFile = fullfile(output_path, '.last_output_path.txt');
            fid = fopen(lastOutputFile, 'w');
            if fid ~= -1
                fprintf(fid, '%s', output_path);
                fclose(fid);
            end
            
            correct_patients = correct_patients + 1;
            disp(['Processing completed successfully for patient: ' current_patient_id]);
        catch ME
            disp(['Error processing patient: ' current_patient_id]);
            disp(ME.message);
            if ~isempty(ME.stack)
                for k = 1:length(ME.stack)
                    fprintf('  File: %s, Line: %d, Function: %s\n', ...
                        ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
                end
            end
            failed_patients = failed_patients + 1;
        end
    end
    
    disp(['Total patients: ' num2str(total_patients)]);
    disp(['Correct patients: ' num2str(correct_patients)]);
    disp(['Failed patients: ' num2str(failed_patients)]);
end

% Nested helper function to process a single patient
function processSinglePatient(path_to_data, eICAB_path, output_path, use_eicab_whole_brain)
    % Load data
    [data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);
    
    % Perform label transfer and preprocessing
    [correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, path_to_data, data_struct, use_eicab_whole_brain);
    
    % Generate LOCs
    [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);
    
    % Compute time-resolved Complex Difference cross-sections
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
    
    % Save data for qvt+
    generateQVTplus(correspondenceDict, LOCs, output_path);
end
end