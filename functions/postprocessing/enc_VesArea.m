function [Area] = enc_VesArea(PI_scat,data_struct,Labels,path2data,params)
    % Load Necessities from labels and processed data
    load(fullfile(path2data,'RawPITC.mat'));
    %%
    BranchList=data_struct.branchList;
    BranchList=[BranchList [1:length(BranchList)]'];
    FlowData=data_struct.flowPulsatile_val;
    roots=[1,2,9]; %L_ICA R_ICA BA
    CoWA=[3,5,4,6,7,8];% LMCA LACA RMCA RACA LPCA RPCA
    Area=zeros([1 9]);
    for i=1:3
        LOC = Labels{roots(i),3};
        LOC=str2num(LOC);
        if length(LOC)>=1
            ves = Labels{roots(i),2};
            ves = str2num(ves);
            if length(ves)>1
                ves = Labels{roots(i),3};
                ves = str2num(ves);
                ves = ves(1);
            end
            if length(LOC)>1
                temp=LOC;
                LOC=temp(2);
                ves=temp(1);
            end
            [idx1,~]=find(BranchList(:,4)==ves);
            Data=BranchList(idx1,:);
            [idx2,~]=find(Data(:,5)==LOC);
            VesLoc=Data(idx2,:);
            [idx3,~]=find(PI_scat(:,3,i)==VesLoc(6));
            Area(1,roots(i))=mean(PI_scat((idx3-2):(idx3+2),5,i));
    end
    loci=[3 5 4 6 7 8];
    for i=1:6
        if i>4
            j=3;
        elseif i>2
            j=2;
        else
            j=1;
        end
        ves = Labels{loci(i),2};
        ves = str2num(ves);
        if length(ves)>1
            ves = Labels{roots(i),3};
            ves = str2num(ves);
            ves=ves(1);
        end
        LOC = Labels{loci(i),3};
        LOC=str2num(LOC);
        if length(LOC)>1
            temp=LOC;
            LOC=temp(2);
            ves=temp(1);
        end
        if ves>0
            [idx1,~]=find(BranchList(:,4)==ves);
            Data=BranchList(idx1,:);
            [idx2,~]=find(Data(:,5)==LOC);
            VesLoc=Data(idx2,:);
            [idx3,~]=find(PI_scat(:,3,j)==VesLoc(6));
            Area(1,loci(i))=mean(PI_scat((idx3-2):(idx3+2),5,j));
        end
    end
    if params.SaveData==1
        save(fullfile(path2data,'Area.mat'),'Area')
    end
end