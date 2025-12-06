function [data_struct, imageData] = saveQVTData(directory,area_val,diam_val,branchList,flowPerHeartCycle_val,maxVel_val,velMean_val,nframes,matrix,res,timeres,...
    VENC,segment,PI_val,RI_val,flowPulsatile_val,r, timeMIPcrossection ,MAGcrossection,segmentFull,autoFlow,vTimeFrameave,...
    Planes,bnumMeanFlow,bnumStdvFlow,StdvFromMean,pixelSpace,VplanesAllx,VplanesAlly,VplanesAllz,imageData,caseFilePath,VoxDims,PIvel_val)

    data_struct = [];
    data_struct.directory = directory;
    data_struct.area_val = area_val;
    data_struct.diam_val = diam_val;
    data_struct.branchList = branchList;
    data_struct.flowPerHeartCycle_val = flowPerHeartCycle_val;
    data_struct.maxVel_val = maxVel_val;
    data_struct.velMean_val = velMean_val;
    data_struct.nframes = nframes;
    data_struct.matrix = matrix;
    data_struct.res = res;
    data_struct.timeres = timeres;
    data_struct.VENC = VENC;
    data_struct.segment = segment;
    data_struct.PI_val = PI_val;
    data_struct.RI_val = RI_val;
    data_struct.flowPulsatile_val = flowPulsatile_val;
    data_struct.r = r;
    data_struct.timeMIPcrossection = timeMIPcrossection;
    data_struct.MAGcrossection = MAGcrossection;
    data_struct.segmentFull = segmentFull;
    data_struct.autoFlow = autoFlow; %SD
    data_struct.vTimeFrameave = vTimeFrameave;
    data_struct.Planes = Planes;
    data_struct.bnumMeanFlow = bnumMeanFlow;
    data_struct.bnumStdvFlow = bnumStdvFlow;
    data_struct.StdvFromMean = StdvFromMean;
    data_struct.pixelSpace = pixelSpace;
    data_struct.VoxDims = VoxDims;
    data_struct.PIvel_val =PIvel_val;
    % Time-resolved CD cross-sections (optional, computed separately)
    if isfield(data_struct, 'timeMIPcrossectionTR')
        % Already set, keep it
    else
        data_struct.timeMIPcrossectionTR = []; % Placeholder if not computed
    end
    % Store original affine matrix to preserve image orientation
    if isfield(imageData, 'OriginalAffine')
        data_struct.OriginalAffine = imageData.OriginalAffine;
    else
        warning('OriginalAffine not found in imageData. Creating default affine from VoxDims.');
        % Fallback: create affine from voxel dimensions (but this should not happen)
        data_struct.OriginalAffine = eye(4);
        data_struct.OriginalAffine(1,1) = VoxDims(1);
        data_struct.OriginalAffine(2,2) = VoxDims(2);
        data_struct.OriginalAffine(3,3) = VoxDims(3);
        new_origin = (data_struct.OriginalAffine.dim(1:3) + 1) / 2; % Center of the image
        data_struct.OriginalAffine(1:3, 4) = data_struct.OriginalAffine(1:3, 1:3) * -new_origin';
    end
    
    Vel_Time_Res.VplanesAllx = VplanesAllx; %TR vel planes (uninterped)
    Vel_Time_Res.VplanesAlly = VplanesAlly;
    Vel_Time_Res.VplanesAllz = VplanesAllz;

    save(caseFilePath,'data_struct','Vel_Time_Res','imageData')