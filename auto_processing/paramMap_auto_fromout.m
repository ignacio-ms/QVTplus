function paramMap_auto_fromout(outputDir, varargin)
% paramMap_auto_fromout - Rerun selected steps of QVT+ pipeline from output directory
%
% This script allows rerunning specific steps of the QVT+ pipeline using
% previously generated output files, avoiding time-consuming operations like
% binary segmentation, registration, and label transfer.
%
% Usage:
%   paramMap_auto_fromout(outputDir)  % Rerun all updateable steps
%   paramMap_auto_fromout(outputDir, 'steps', {'LOCs', 'saveVesselData', 'generateQVTplus'})
%   paramMap_auto_fromout(outputDir, 'steps', {'LOCs'})  % Only regenerate LOCs
%   paramMap_auto_fromout(outputDir, 'steps', {'LOCs', 'timeResolvedCD'})
%
% Inputs:
%   outputDir - Path to QVT+ output directory (folder containing
%               qvtData_ISOfix_*.mat, SummaryParamTool.xls,
%               multilabel_QVTseg.nii, etc.)
%
% Optional Name-Value Pairs:
%   'steps' - Cell array of steps to rerun. Valid steps:
%             - 'LOCs': Regenerate LOCs using new algorithms
%             - 'saveVesselData': Update vessel-specific Excel tables and images
%             - 'generateQVTplus': Update LabelsQVT.csv file
%             - 'timeResolvedCD': Recompute time-resolved Complex Difference cross-sections & Registration
%             Default: {'LOCs', 'saveVesselData', 'generateQVTplus'}
%
% Examples:
%   % Rerun all updateable steps
%   paramMap_auto_fromout('/results/sub-001')
%
%   % Only regenerate LOCs (e.g., with new circularity-based selection)
%   paramMap_auto_fromout('/results/sub-001', 'steps', {'LOCs'})
%
%   % Regenerate LOCs and update all tables
%   paramMap_auto_fromout('/results/sub-001', 'steps', {'LOCs', 'saveVesselData', 'generateQVTplus'})

    % Parse input arguments
    if nargin < 1 || isempty(outputDir)
        error('paramMap_auto_fromout:MissingOutputDir', ...
              'outputDir is required. Use --help for usage information.');
    end
    
    if nargin == 2 && (strcmp(varargin{1}, '--help') || strcmp(varargin{1}, '-h') || strcmp(varargin{1}, '-?'))
        help('paramMap_auto_fromout');
        return;
    end
    
    % Parse optional name-value pairs
    p = inputParser;
    validSteps = {'LOCs', 'saveVesselData', 'generateQVTplus', 'timeResolvedCD'};
    defaultSteps = {'LOCs', 'saveVesselData', 'generateQVTplus'};
    
    addRequired(p, 'outputDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'steps', defaultSteps, @(x) iscell(x) && all(ismember(x, validSteps)));
    
    parse(p, outputDir, varargin{:});
    
    outputDir = char(p.Results.outputDir);
    stepsToRun = p.Results.steps;
    
    % Validate output directory
    if ~exist(outputDir, 'dir')
        error('paramMap_auto_fromout:DirectoryNotFound', ...
              'Output directory not found: %s', outputDir);
    end
    
    % Check for required files
    matInfo = dir(fullfile(outputDir, 'qvtData_ISOfix_*.mat'));
    if isempty(matInfo)
        error('paramMap_auto_fromout:MissingData', ...
              'No qvtData_ISOfix_*.mat file found in %s.', outputDir);
    end
    
    multilabelPath = fullfile(outputDir, 'multilabel_QVTseg.nii');
    if ~exist(multilabelPath, 'file')
        warning('paramMap_auto_fromout:MissingMultilabel', ...
                'multilabel_QVTseg.nii not found. Some steps may fail.');
    end
    
    % Display what we're doing
    disp(['================================']);
    disp(['paramMap_auto_fromout']);
    disp(['Output directory: ' outputDir]);
    disp(['Steps to run: ' strjoin(stepsToRun, ', ')]);
    disp(['================================']);
    
    % Load existing data
    [~, newestIdx] = max([matInfo.datenum]);
    matFile = fullfile(outputDir, matInfo(newestIdx).name);
    
    disp('Loading existing data...');
    dataVars = load(matFile, 'data_struct', 'Vel_Time_Res', 'imageData');
    data_struct = dataVars.data_struct;
    Vel_Time_Res = dataVars.Vel_Time_Res;
    
    if isfield(dataVars, 'imageData')
        imageData = dataVars.imageData;
    else
        imageData = [];
        if ismember('timeResolvedCD', stepsToRun)
            warning('paramMap_auto_fromout:MissingImageData', ...
                    'imageData not found. Cannot compute time-resolved CD.');
            stepsToRun = setdiff(stepsToRun, {'timeResolvedCD'});
        end
    end
    
    % Rebuild correspondence dictionary from existing files
    try
        disp('Rebuilding correspondence dictionary...');
        [correspondenceDict, multiQVT] = generateCorrespondenceDict(outputDir, data_struct);
    catch ME
        error('paramMap_auto_fromout:CorrespondenceFailed', ...
              'Unable to rebuild correspondence dictionary: %s', ME.message);
    end
    
    % Run selected steps
    LOCs = [];
    
    if ismember('LOCs', stepsToRun)
        disp('--------------------------------');
        disp('Step: Regenerating LOCs...');
        try
            [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);
            disp('LOCs regenerated successfully.');
        catch ME
            warning('paramMap_auto_fromout:LOCsFailed', ...
                    'Failed to regenerate LOCs: %s', ME.message);
            if ~isempty(ME.stack)
                for k = 1:length(ME.stack)
                    fprintf('  File: %s, Line: %d, Function: %s\n', ...
                        ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
                end
            end
        end
    end
    
    if ismember('timeResolvedCD', stepsToRun) && ~isempty(imageData)
        disp('--------------------------------');
        disp('Step: Computing time-resolved Complex Difference cross-sections...');
        try
            % Determine path_to_data from output directory structure
            % Try to infer from output path (assumes standard structure)
            path_to_data = [];
            if exist(fullfile(outputDir, '..', '..', 'DATA', 'Nifti'), 'dir')
                % Batch mode structure
                [~, patient_id] = fileparts(outputDir);
                path_to_data = fullfile(outputDir, '..', '..', 'DATA', 'Nifti', patient_id, '4DFlow');
                if ~exist(path_to_data, 'dir')
                    path_to_data = [];
                end
            end
            
            if isempty(path_to_data)
                % Prompt user or use a default - for now, we'll try to use imageData
                path_to_data = outputDir; % Fallback
                warning('paramMap_auto_fromout:InferringPath', ...
                        'Could not infer path_to_data. Using output directory as fallback.');
            end
            
            timeMIPcrossectionTR = computeTimeResolvedCD(data_struct, imageData, path_to_data, outputDir);
            % Add to data_struct for saving
            data_struct.timeMIPcrossectionTR = timeMIPcrossectionTR;
            
            % Update the saved .mat file
            savedData = load(matFile);
            savedData.data_struct.timeMIPcrossectionTR = timeMIPcrossectionTR;
            save(matFile, '-struct', 'savedData', '-v7.3');
            disp('Time-resolved CD data saved successfully.');
        catch ME
            warning('paramMap_auto_fromout:TimeResolvedCDFailed', ...
                    'Failed to compute time-resolved CD: %s', ME.message);
            if ~isempty(ME.stack)
                for k = 1:length(ME.stack)
                    fprintf('  File: %s, Line: %d, Function: %s\n', ...
                        ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
                end
            end
        end
    end
    
    if ismember('saveVesselData', stepsToRun)
        disp('--------------------------------');
        disp('Step: Updating vessel-specific data...');
        try
            if isempty(LOCs)
                % Load existing LOCs if not regenerated
                % Try to load from generateCorrespondenceDict (which rebuilds from files)
                [~, tempLOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);
                LOCs = tempLOCs;
            end
            
            saveVesselData(LOCs, data_struct, outputDir);
            disp('Vessel-specific data updated successfully.');
        catch ME
            warning('paramMap_auto_fromout:SaveVesselDataFailed', ...
                    'Failed to update vessel-specific data: %s', ME.message);
            if ~isempty(ME.stack)
                for k = 1:length(ME.stack)
                    fprintf('  File: %s, Line: %d, Function: %s\n', ...
                        ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
                end
            end
        end
    end
    
    if ismember('generateQVTplus', stepsToRun)
        disp('--------------------------------');
        disp('Step: Updating QVT+ labels file...');
        try
            if isempty(LOCs)
                % Load existing LOCs if not regenerated
                [~, tempLOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);
                LOCs = tempLOCs;
            end
            
            generateQVTplus(correspondenceDict, LOCs, outputDir);
            disp('QVT+ labels file updated successfully.');
        catch ME
            warning('paramMap_auto_fromout:GenerateQVTplusFailed', ...
                    'Failed to update QVT+ labels: %s', ME.message);
            if ~isempty(ME.stack)
                for k = 1:length(ME.stack)
                    fprintf('  File: %s, Line: %d, Function: %s\n', ...
                        ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
                end
            end
        end
    end
    
    disp(['================================']);
    disp(['Processing completed.']);
    disp(['================================']);
end
