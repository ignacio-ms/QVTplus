function generateQVTplus(correspondenceDict, LOCs, output_path)
    % Function to generate a table of artery labels and LOCs, and save it as a CSV file.
    %
    % Parameters:
    %   correspondenceDict - Structure containing correspondence information
    %   LOCs - Structure containing LOC information
    %   outputFileName - Name of the output CSV file

    % Define artery key mapping and names
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

    % Initialize cell arrays for table columns
    arteryList = {};
    labelList = {};
    locList = {};

    % Loop through each artery name and extract data
    for i = 1:numel(arteryNames)
        artery = arteryNames{i};
        internalKey = arteryKeyMap.(strrep(artery, ' ', '_'));

        % Initialize strings for label and LOC
        labelStr = ' ';
        locStr = ' ';

        % Retrieve label data if available
        if ~isempty(internalKey) && isfield(correspondenceDict, internalKey)
            labelData = correspondenceDict.(internalKey);
            labelStr = sprintf('[%s]', strjoin(string(labelData), ', '));
        end

        % Retrieve LOC data if available
        if ~isempty(internalKey) && isfield(LOCs, internalKey)
            locData = LOCs.(internalKey);
            if numel(locData) > 1
                locStr = sprintf('[%d, %d]', locData(1), locData(2));
            else
                locStr = sprintf('%d', locData(1));
            end
        end

        % Append to the lists
        arteryList{end+1} = artery; %#ok<AGROW>
        labelList{end+1} = labelStr; %#ok<AGROW>
        locList{end+1} = locStr; %#ok<AGROW>
    end

    % Create the table
    T = table(arteryList', labelList', locList', 'VariableNames', {'Artery', 'Label', 'Loc'});

    % Write the table to the output CSV file
    writetable(T, fullfile(output_path, 'LabelsQVT.csv'));
end
