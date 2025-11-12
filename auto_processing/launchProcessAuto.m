function launchProcessAuto(path_to_data, eicab_path, output_path)
    % Load data
    [data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

    % Perform label transfer and preprocessing
    [correspondenceDict, multiQVT] = performLabelTransfer(eicab_path, output_path, imageData, path_to_data, data_struct);

    % Generate LOCs
    [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

    % Save vessel-specific data automatically
    saveVesselData(LOCs, data_struct, output_path);

    %Save data for qvt+
    generateQVTplus(correspondenceDict, LOCs, output_path)

    disp('Processing completed successfully.');
end
