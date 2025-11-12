function varargout = correspondence_funcs(funcName, varargin)
    % Dispatcher function
    switch funcName
        case 'buildLabelOccurrences'
            varargout{1} = buildLabelOccurrences(varargin{:});
        case 'resolveLabelMappings'
            varargout{1} = resolveLabelMappings(varargin{:});
        case 'removeDuplicateEntries'
            varargout{1} = removeDuplicateEntries(varargin{:});
        otherwise
            error('Unknown function name: %s', funcName);
    end
end


function labelOccurrences = buildLabelOccurrences(correspondenceDict)
    labelOccurrences = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    keys = fieldnames(correspondenceDict);

    for k = 1:numel(keys)
        key = keys{k};
        labelList = correspondenceDict.(key);

        for label = labelList'
            if ~isKey(labelOccurrences, label)
                labelOccurrences(label) = struct('keys', {{}}, 'counts', []);
            end

            entry = labelOccurrences(label);
            if ~ismember(key, entry.keys)
                entry.keys{end+1} = key;
                entry.counts(end+1) = 0;
            end

            keyIdx = find(strcmp(entry.keys, key));
            entry.counts(keyIdx) = entry.counts(keyIdx) + 1;

            labelOccurrences(label) = entry;
        end
    end
end

function correspondenceDict = resolveLabelMappings(correspondenceDict, labelOccurrences)
    % Iterate through keys of the labelOccurrences map
    labelKeys = keys(labelOccurrences); % Get cell array of keys
    for i = 1:numel(labelKeys)
        label = labelKeys{i}; % Extract the current key
        entry = labelOccurrences(label);

        % Ignore 'good_lab_100'
        ignoreIndex = strcmp(entry.keys, 'good_lab_100');
        validCounts = entry.counts(~ignoreIndex);
        validKeys = entry.keys(~ignoreIndex);

        if isempty(validCounts)
            continue;
        end

        % Find the key with the highest count
        [~, maxIdx] = max(validCounts);
        keepKey = validKeys{maxIdx};

        % Remove the label from all other keys
        for k = 1:numel(entry.keys)
            key = entry.keys{k};
            if ~strcmp(key, keepKey)
                correspondenceDict.(key) = correspondenceDict.(key)(correspondenceDict.(key) ~= label);
            end
        end
    end
end

function correspondenceDict = removeDuplicateEntries(correspondenceDict)
    keys = fieldnames(correspondenceDict);
    for k = 1:numel(keys)
        key = keys{k};
        correspondenceDict.(key) = unique(correspondenceDict.(key));
    end
end