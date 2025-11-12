function [data_struct, imageData] = loadPreprocessedData(path_to_data, output_path)
    matFiles = dir(fullfile(path_to_data, '*.mat'));
    if ~isempty(matFiles)
        disp('Loading preprocessed files...');
        data_struct = load(fullfile(path_to_data, matFiles(1).name));
    else
        disp('Loading raw data...');

        % Check if the data is in DICOM or NIFTI format
        
        % if there are *.nii or *.nii.gz files in the directory, use the loadNII_auto function
        niiFiles = dir(fullfile(path_to_data, '**/*.nii*'));

        if ~isempty(niiFiles)
            [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
            maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
            VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
            timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
            bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullEx,autoFlow,pixelSpace, VoxDims, PIvel_val] = loadNII_auto(path_to_data);
        elseif exist([path_to_data filesep],'dir')
            [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
            maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
            VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
            timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
            bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullEx,autoFlow,pixelSpace, VoxDims, PIvel_val] = loadDCM(directory);
        else
            error('No valid data found in the specified path.');
        end

        % create the output directory if it does not exist
        if ~exist(output_path, 'dir')
            mkdir(output_path);
        end

        % Save the processed data
        time = datestr(now, 'ddmmmyyyy_HHMM');
        caseFilePath = fullfile(output_path, ['qvtData_ISOfix_' time '.mat']);
        [data_struct, imageData] = saveQVTData(caseFilePath,area_val,diam_val,branchList,flowPerHeartCycle_val,maxVel_val,velMean_val,nframes,matrix,res,timeres,...
        VENC,segment,PI_val,RI_val,flowPulsatile_val,r, timeMIPcrossection ,MAGcrossection,segmentFull,autoFlow,vTimeFrameave,...
        Planes,bnumMeanFlow,bnumStdvFlow,StdvFromMean,pixelSpace,VplanesAllx,VplanesAlly,VplanesAllz,imageData,caseFilePath,VoxDims,PIvel_val);
         
    end
end
