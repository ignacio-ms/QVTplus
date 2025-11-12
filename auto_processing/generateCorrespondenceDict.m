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