% ProcessPITC loads the QVT saved datastruct and computes the pulsatility
% as a function of length. For fitting PITC, see FitPITC.
%
% All that's needed is the path to QVT saved data which must include
% the LabelsQVT.csv
% You must also initialise a search distance for the connctivity algorithm.
% The returned connectivity should be checked with the basilar root for any
% missed PCA's, and increased if necessary.
%
% Outputs: saved the processed connectivity maps, and pulsatility and
% length matrices.
%
% Used by: None, BUT, after use, run fitPITC. The functions are separated
% for no reason haha.
%
% Dependencies: compute_conn.m, compute_length.m
%% Initialise
clear;clc;
subject='sub-014';
path2data=strcat('C:\Users\sdem348\Desktop\Dempsey2023MultiRes_Cohort\derivatives\QVT\',subject); % just to the root folder, can only need this if not using bids
SearchDist=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Don't change below %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Or do
%% Load Necessities from labels and processed data
DataName = dir(fullfile(path2data,'*.mat'));
for i=1:length(DataName)
    if length(DataName(i).name) >6
        if strcmp(DataName(i).name(1:7),'qvtData')
            LOC=i;
            break
        end
    end
end
load(fullfile(path2data,DataName(LOC).name));
clear LOC DataName
Labels=readLabels(path2data);
%% CONNECTIVITY SECTION: RUN THIS ITERATIVELY IF THE MATS CONNECTIVITIES DON'T MATCH BASILAR ROOT
ICA_l = Labels{1,2} ;ICA_l = str2num(ICA_l);
ICA_r = Labels{2,2} ;ICA_r = str2num(ICA_r);
BA = Labels{9,2}    ;BA = str2num(BA);
Exclude=Labels{10,2};Exclude = str2num(Exclude);
LR=Labels{10,3}     ;LR = str2num(LR);
startroots = {ICA_l,ICA_r,BA};
%startroots = {ICA_r,ICA_l,BA}; %SWITCHED LR SIDE
BranchList=data_struct.branchList;
BranchList=[BranchList [1:length(BranchList)]'];
for i=1:length(Exclude)
    [idx,~]=find(BranchList(:,4)==Exclude(i));
    BranchList(idx,:)=[];
end
% Compute Connectivity below
Mats=compute_conn(startroots,BranchList,SearchDist,LR);
save(fullfile(path2data,'ConnectivityMaps.mat'),'Mats');
%% LENGTH SECTION
load(fullfile(path2data,'ConnectivityMaps.mat'));
BranchList=data_struct.branchList;
BranchList=[BranchList [1:length(BranchList)]'];
PI=data_struct.PI_val;
Quality=data_struct.StdvFromMean;
MaxVel=data_struct.PIvel_val;
VoxDims=data_struct.VoxDims;
% Apply the voxel dimensions to the matrix location of the branch lists
BranchList(:,1)=VoxDims(1).*BranchList(:,1);
BranchList(:,2)=VoxDims(2).*BranchList(:,2);
BranchList(:,3)=VoxDims(3).*BranchList(:,3);
% Process Velocity Based Pulsatility (default output is flow based)
PIvel_val=zeros([length(MaxVel(:,1)) 1]);
for i=1:length(MaxVel(:,1))
    PIvel_val(i)=(max(MaxVel(i,:))-min(MaxVel(i,:)))./(mean(MaxVel(i,:)));
end
% Compute Lengths Below
[PI_scat,PIvel_scat,G] = compute_length(Mats,BranchList,PI,PIvel_val,Quality);
%save(fullfile(path2data,'RawDamping.mat'),'PI_scat',"PIvel_scat",'G');