%% Running a Directory
%==================================================
% Running a BIDS processed derivates folder for PI, PI_DF, PI_TC, PWV is below
%==================================================
% This function is designed to run assuming you store your data in BIDS 
% format. If not, modify to run a folder or what not, BIDS is much easier
clear;clc;
path2bids='C:\Users\u149879\Desktop\trial\dicom2';
params=init_params(path2bids);
params.SaveData=1;%Save Interim Data
params.PltFlag=1; %Make plots
params.rerun=1; %esto es solo por si queremos hacer el linkage analysis por cada caso
Results=struct;
for i=1:2
    subject=strcat('sub-',num2str(i,'%03.0f'));
    params.subject=subject;
    path2data=path2bids;
    fprintf(strcat('Now Processing Case:',num2str(i),'\n'))
    if ~exist(string(fullfile(path2data,subject)), 'dir')
        mkdir(string(fullfile(path2data,subject)))
    end
    %string(fullfile(params.data_dir,'derivatives\QVT',subject));
    if exist(path2data,'dir')==7
        [data_struct,Labels,~]=import_QVTplusData(path2data);
        [PI_scat]=enc_PITC_process(data_struct,Labels,path2data,params,[]);
        [PITC,globPI] = enc_PITC_fit(PI_scat,path2data,params);
        [DF,PI,LocFlows,FlowErr,D] = enc_PI_DF_Flows(PI_scat,data_struct,Labels,path2data,params);
        [Area] = enc_VesArea(PI_scat,data_struct,Labels,path2data,params);
        [time,Flow,FlowErrVes]=enc_HQVesselFlows(data_struct,Labels,params);
        [PWV,R]=enc_PWV(data_struct,PI_scat,time,Labels,params);
        Results.PITC(i,:)=PITC(1,:);
        Results.PI(i,:)=[globPI PITC(end,:) PI];
        Results.DF(i,:)=DF;
        Results.Flow(i,:)=[mean(LocFlows)];
        Results.PWV(i,:)=PWV;
        Results.PWVr(i,:)=R;
        Results.Area(i,:)=Area;
    else
    end
end
if ~exist(string(fullfile(path2data,'population')), 'dir')
        mkdir(string(fullfile(path2data,'population')))
end
save(fullfile(params.data_dir,'population\ResultsMR.mat'),"Results")

% %% Running 
% %==================================================
% % Running a single instance of the post processing for PI, PI_DF, PI_TC, PWV is below
% %==================================================
% clear;clc;
% path2data='C:\Users\sdem348\Documents\MATLAB\CHM\QVTplus\testdata'; %should have the  LabelsQVT and qvtData inside
% params=init_params(path2data);
% params.PltFlag=1; %Plot your results
% params.SaveData=1; %Save results and plots
% params.rerun=1; %rerun PITCH algorithm, or use stored data
% Results=struct;
% [data_struct,Labels]=import_QVTplusData(path2data);
% [PI_scat]=enc_PITC_process(data_struct,Labels,path2data,params,[]);
% [PITC,globPI] = enc_PITC_fit(PI_scat,path2data,params);
% [DF,PI,LocFlows,LocFlowErr] = enc_PI_DF_Flows(PI_scat,data_struct,Labels,path2data,params);
% [time,Flow,FlowErr]=enc_HQVesselFlows(data_struct,Labels,params);
% close all
% Results.PITC=PITC(1,:);
% Results.PI=[globPI PITC(end,:) PI];
% Results.DF=DF;
% Results.Flow=[mean(LocFlows)];
% save(fullfile(path2data,'ResultsALL.mat'),"Results")
