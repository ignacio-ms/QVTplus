function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullEx,autoFlow,pixelSpace, VoxDims, PIvel_val] = loadNII_auto(directory)
% loadNII loads NIFTI files and then processes them.
%
% It returns too much to discuss, but it basically passes through all the
% processed data to paramMap.
%
% Outputs: Everything?
%
% Used by: autoCollectFlow.m, and any separate functions to compute any
% saved data (PITC codes, Damping codes etc)
%
% Notes: The DICOM data may be rotated 180 (reverse flows), this may be
% allieviated by rotating the matrices using something like 
% %MAG=imrotate(MAG(:,:,:),rotImAngle);
% This should be done for all loaded matrices, HOWEVER, this was only
% required on Siemens data for some reason and may not always be needed.
%
% Dependencies: ShuffleDCM.m and all QVT processing codes
%% Initialization
clc
BGPCdone=0; %0=do backgroun correction, 1=don't do background correction.
%VENC = 700; %may change depending on participant
autoFlow=1; %if you want automatically extracted BC's and flow profiles 0 if not.

filetype = 'nii';

% Check if 'scans' subfolder exists
if isfolder(fullfile(directory, 'scans'))
    base_dir = fullfile(directory, 'scans');
    is_scans_present = true;
else
    base_dir = directory;
    is_scans_present = false;
end

% Find direction folders with different approaches depending on structure
if is_scans_present
    % Original approach with scans subfolder
    folders = dir(base_dir);
    folders(ismember({folders.name}, {'.', '..'})) = [];
    folders = {folders([folders.isdir]).name};
    
    folderAP = folders(~cellfun('isempty', regexp(folders, 'AP')));
    folderAP = folderAP{1};
    folderRL = folders(~cellfun('isempty', regexp(folders, 'RL')));
    folderRL = folderRL{1};
    folderFH = folders(~cellfun('isempty', regexp(folders, 'FH')));
    folderFH = folderFH{1};
else
    % Directly search for direction folders in the main directory
    all_dirs = dir(base_dir);
    all_dirs = all_dirs([all_dirs.isdir]);
    dir_names = {all_dirs.name};
    
    folderAP = dir_names(~cellfun('isempty', regexp(dir_names, 'AP')));
    folderAP = folderAP{1};
    folderRL = dir_names(~cellfun('isempty', regexp(dir_names, 'RL')));
    folderRL = folderRL{1};
    folderFH = dir_names(~cellfun('isempty', regexp(dir_names, 'FH')));
    folderFH = folderFH{1};
end

% Get paths for each direction
ap_path = fullfile(base_dir, folderAP);
rl_path = fullfile(base_dir, folderRL);
fh_path = fullfile(base_dir, folderFH);

% Check if NIFTI subfolder exists in each direction folder
if isfolder(fullfile(ap_path, 'NIFTI'))
    ap_path = fullfile(ap_path, 'NIFTI');
    rl_path = fullfile(rl_path, 'NIFTI');
    fh_path = fullfile(fh_path, 'NIFTI');
end

% Load magnitude JSON file
json_files = dir(fullfile(ap_path, '*.json'));
json_mag = fileread(fullfile(json_files(1).folder, json_files(1).name));
json_mag = jsondecode(json_mag);

% Load magnitude volume
magvol = dir(fullfile(ap_path, '*.nii.gz'));
magvol = spm_vol(fullfile(magvol(1).folder, magvol(1).name));
mag = flip(spm_read_vols(magvol), 3);

% Load phase volumes for each direction
vxvol = dir(fullfile(ap_path, '*_ph.nii.gz'));
vxvol = spm_vol(fullfile(vxvol(1).folder, vxvol(1).name));
vx = flip(spm_read_vols(vxvol), 3);

vyvol = dir(fullfile(rl_path, '*_ph.nii.gz'));
vyvol = spm_vol(fullfile(vyvol(1).folder, vyvol(1).name));
vy = flip(spm_read_vols(vyvol), 3);

vzvol = dir(fullfile(fh_path, '*_ph.nii.gz'));
vzvol = spm_vol(fullfile(vzvol(1).folder, vzvol(1).name));
vz = flip(spm_read_vols(vzvol), 3);

[a,c,b,d] = size(vx);
v = zeros([a,c,b,3,d],'single');

% velocities are in cm/s, convert to mm/s

v(:,:,:,2,:)=-squeeze(vx(:,:,:,:))*10;
v(:,:,:,1,:)=squeeze(vy(:,:,:,:))*10;
v(:,:,:,3,:)=-squeeze(vz(:,:,:,:))*10;

vMean = mean(v,5);

MAG = mean(mag,4);

nframes = length(magvol);

% are these needed?
VENC=700; %VENC is 700 mm/s, 70 cm/s
timeres = 60/(60*nframes);%INFO.NominalInterval/nframes; %temporal resolution (ms) ????

matrix = magvol(1).dim;
VoxDims = sqrt(sum(magvol(1).mat(1:3,1:3).^2));
res = VoxDims(1);
slicespace = VoxDims(3);

%% Import Complex Difference
timeMIP = calc_angio(MAG, vMean, VENC);

%% Manual Background Phase Correction (if necessary)
back = zeros(size(vMean),'single');
if ~BGPCdone
    [poly_fitx,poly_fity,poly_fitz] = unattended_background_phase_correction(MAG,vMean(:,:,:,1),vMean(:,:,:,2),vMean(:,:,:,3),2,0.11,0.06);
    disp('Correcting data with polynomial');
    xrange = single(linspace(-1,1,size(MAG,1)));
    yrange = single(linspace(-1,1,size(MAG,2)));
    zrange = single(linspace(-1,1,size(MAG,3)));
    [Y,X,Z] = meshgrid(yrange,xrange,zrange);
    % Get poly data and correct average velocity for x,y,z dimensions
    back(:,:,:,1) = single(evaluate_poly(X,Y,Z,poly_fitx));
    back(:,:,:,2) = single(evaluate_poly(X,Y,Z,poly_fity));
    back(:,:,:,3) = single(evaluate_poly(X,Y,Z,poly_fitz));
    vMean = vMean - back;
    for f=1:nframes
        v(:,:,:,:,f) = v(:,:,:,:,f) - back;
    end 
    clear X Y Z poly_fitx poly_fity poly_fitz xrange yrange zrange
end
%% Find optimum global threshold for total branch segmentation
step = 0.001; %step size for sliding threshold
UPthresh = 0.8; %max upper threshold when creating Sval curvature plot
SMf = 10;
shiftHM_flag = 1; %flag to shift max curvature by FWHM
medFilt_flag = 1; %flag for median filtering of CD image
own_seg = 0;
if own_seg == 0
    [~,segment] = slidingThreshold(timeMIP,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
    areaThresh = round(sum(segment(:)).*0.005); %minimum area to keep
    conn = 6; %connectivity (i.e. 6-pt)
    segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes
else
    NII = spm_vol('/home/afernandezpe/borrar/4DFlow_reg/tests_QVT/rBMRI198894_4DQflowNeuro_AP_1803_tAvgProj_eICAB_WB.nii.gz');
    segment = spm_read_vols(NII);
    segment(isnan(segment)) = 0;
end
% save raw (cropped) images to imageData structure (for Visual Tool)
imageData.MAG = MAG;
imageData.CD = timeMIP; 
imageData.V = vMean;
imageData.Segmented = segment;
imageData.Header = json_mag;

%save("C:\Users\u149879\Desktop\trial\original_ABI\ABI_nifti",imageData,'-mat')

%% Feature Extraction
% Get trim and create the centerline data
sortingCriteria = 3; %sorts branches by junctions/intersects 
spurLength = 8; %minimum branch length (removes short spurs)
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment);

%% You can load another segmentation here if you want which will overlap on images
Exseg=segment; %for now, dummy copy, can do feature extraction
%[~,~,branchList2,~] = feature_extraction(sortingCriteria,spurLength,vMean,logical(JSseg),[]);

%% SEND FOR PROCESSING
% Flow parameter calculation, bulk of code is in paramMap_parameters.m
SEG_TYPE = 'thresh'; %kmeans or thresh
if strcmp(SEG_TYPE,'kmeans')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_kmeans(filetype,branchList,matrix,timeMIP,vMean, ...
    back,BGPCdone,directory,nframes,res,MAG, v,slicespace);
elseif strcmp(SEG_TYPE,'thresh')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes,pixelSpace,...
        segmentFullEx,PIvel_val] ...
        = paramMap_params_threshS(filetype,branchList,matrix,timeMIP,vMean, ...
    back,BGPCdone,directory,nframes,res,MAG, v,slicespace,Exseg);
else
    disp("Incorrect segmentation type selected, please select 'kmeans' or 'thresh'");
end 
return