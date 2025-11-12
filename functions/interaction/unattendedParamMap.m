function unattendedParamMap(directory)
% Unattended version of paramMap that takes inputs rather than using GUI
% Inputs:
%   directory - Path containing 4D flow data
%   

versionNum = 'v1-2'; %paramMap Version

% Initialize global variables as in original
global branchList Planes hfull p branchLabeled Ntxt nframes res matrix VENC
global AveAreaBranch LogPoints fullCData area_val flowPerHeartCycle_val
global PI_val diam_val maxVel_val RI_val flowPulsatile_val timeres segment
global r timeMIPcrossection segmentFull vTimeFrameave velMean_val
global MAGcrossection bnumMeanFlow bnumStdvFlow StdvFromMean
global VplanesAllx VplanesAlly VplanesAllz imageData
global segmentFullJS autoFlow pixelSpace VoxDims PIvel_val

% Load data based on directory contents
if exist([directory filesep 'scans'],'dir') 
    [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val,...
        maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val,...
        VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r,...
        timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection,imageData,...
        bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullJS,autoFlow,pixelSpace,...
        VoxDims,PIvel_val] = loadNII(directory);
elseif exist([directory filesep],'dir')
    [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val,...
        maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val,...
        VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r,...
        timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection,imageData,...
        bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullJS,autoFlow,pixelSpace,...
        VoxDims,PIvel_val] = loadDCM(directory);
end

time = datestr(now);
saveState = [time(1:2) time(4:6) time(10:11) '_' time(13:14) time(16:17) '_' versionNum];
caseFilePath = fullfile(directory,['qvtData_ISOfix_' saveState '.mat']);

% save the data
saveQVTData(directory,area_val,diam_val,branchList,flowPerHeartCycle_val,maxVel_val,velMean_val,nframes,matrix,res,timeres,...
VENC,segment,PI_val,RI_val,flowPulsatile_val,r, timeMIPcrossection ,MAGcrossection,segmentFull,segmentFullJS,autoFlow,vTimeFrameave,...
Planes,bnumMeanFlow,bnumStdvFlow,StdvFromMean,pixelSpace,VplanesAllx,VplanesAlly,VplanesAllz,imageData,caseFilePath,VoxDims,PIvel_val)

% This will be the name used for the Excel file
finalFolder = regexp(directory,filesep,'split');
SummaryName = [finalFolder{end} '_qvtData_ISOfix_' saveState];
warning off
mkdir(directory,SummaryName); %makes directory if it already exists

% Where to save data images and excel summary files
SavePath = [directory filesep SummaryName];
if autoFlow == 1
    autoCollectFlow(directory)
    %autoCollectDamping(directory)
end


% Create excel files save summary data
vesselNames = {'Left ICA Cavernous';'Left ICA Cervical';'Right ICA Cavernous'; ...
    'Right ICA Cervical';'Basilar';'Left MCA';'Right MCA';'Left PCA';...
    'Right PCA';'SS sinus';'Straight sinus';'Left Transverse';...
    'Right Transverse';'Left VA';'Right VA';'Left ACA';'Right ACA'; ...
    'Misc1'; 'Misc2'};

col_header = ({'Vessel Label', 'Centerline Point', 'Notes',['Max Velocity < ' num2str(VENC) 'cm/s'], ...
    'Mean Flow ml/s','Pulsatility Index','Branch Label'});
xlwrite([SavePath filesep 'SummaryParamTool.xls'],col_header,'Summary_Centerline','A1');
xlwrite([SavePath filesep 'SummaryParamTool.xls'],vesselNames,'Summary_Centerline','A2');


end