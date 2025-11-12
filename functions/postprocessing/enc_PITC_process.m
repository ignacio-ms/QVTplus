function [PI_scat] = enc_PITC_process(data_struct,Labels,path2data,params,SDist)
    %% Load Necessities from labels and processed data
    if ~isempty(path2data)
        if params.rerun == 0 %Don't rerun unless you want to for some reason
            if isfile(fullfile(path2data,'RawPITC.mat'))
                ProcessFlag=0;
                load(fullfile(path2data,'RawPITC.mat'));
            else
                ProcessFlag=1;
            end
        else
            ProcessFlag=1;
        end
    else
        ProcessFlag=1;
    end
    if ProcessFlag==1
        ICA_l = Labels{1,2} ;ICA_l = str2num(ICA_l);
        ICA_r = Labels{2,2} ;ICA_r = str2num(ICA_r);
        BA = Labels{9,2}    ;BA = str2num(BA);
        Exclude=Labels{10,2};Exclude = str2num(Exclude);
        LR=Labels{10,3}     ;LR = str2num(LR);
        startroots = {ICA_l,ICA_r,BA};
        if isempty(SDist)
            ds = DataStorage();  
            my_app_handle=PITCinteractive2(data_struct, startroots, LR, Exclude,ds);
            uiwait(my_app_handle.UIFigure);
            SearchDist=ds.dataArea.SearchDist;
            PI_scat=ds.dataArea.PI_scat;
            G=ds.dataArea.G;
            Mats=ds.dataArea.Mats;
            if params.SaveData==1
                save(fullfile(path2data,'RawPITC.mat'),'PI_scat','G','Mats','SearchDist');
            end
        else
            BranchList=data_struct.branchList;
            BranchList=[BranchList [1:length(BranchList)]'];
            for i=1:length(Exclude)
                [idx,~]=find(BranchList(:,4)==Exclude(i));
                BranchList(idx,:)=[];
            end
            Mats=compute_conn(startroots,BranchList,SDist,LR);
            save(fullfile(path2data,'ConnectivityMaps.mat'),'Mats');
            %% LENGTH SECTION
            load(fullfile(path2data,'ConnectivityMaps.mat'));
            BranchList=data_struct.branchList;
            BranchList=[BranchList [1:length(BranchList)]'];
            PI=data_struct.PI_val;
            Quality=data_struct.StdvFromMean;
            Area=data_struct.area_val;
            VoxDims=data_struct.VoxDims;
            % Apply the voxel dimensions to the matrix location of the branch lists
            BranchList(:,1)=VoxDims(1).*BranchList(:,1);
            BranchList(:,2)=VoxDims(2).*BranchList(:,2);
            BranchList(:,3)=VoxDims(3).*BranchList(:,3);
            % Compute Lengths Below
            [PI_scat,G,~] = compute_length(Mats,BranchList,PI,Area,Quality);
            SearchDist=SDist;
            if params.SaveData==1
                save(fullfile(path2data,'RawPITC.mat'),'PI_scat','G','Mats','SearchDist');
            end
        end
    end
end