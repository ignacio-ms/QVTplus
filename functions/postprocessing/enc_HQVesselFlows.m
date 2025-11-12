function [time,Flow,FlowErr]=enc_HQVesselFlows(data_struct,Labels,params)
    BranchList=data_struct.branchList;
    BranchList=[BranchList [1:length(BranchList)]'];
    Quality=data_struct.StdvFromMean;
    BranchList=[BranchList Quality];
    Flows=data_struct.flowPulsatile_val;
    [~,frames]=size(Flows);
    time=(data_struct.timeres:data_struct.timeres:(data_struct.timeres*(frames)));
    Flow=zeros([frames,9]);
    FlowErr=zeros([frames,9]);
    for ves=1:9
        Vessel=Labels{ves,2};Vessel=str2num(Vessel);
        Data=[];
        if ~isempty(Vessel) %In case standard vessel doesn't exist (ACAs from 1 ICA)
            for vessnum=1:length(Vessel)
                [idx1,~]=find(BranchList(:,4)==Vessel(vessnum));
                Temp=BranchList(idx1,:);
                [idx2,~]=find(Temp(:,7)>=params.thresh);
                Data=[Data;Temp(idx2,:)];
            end
            HQFlows=Flows(Data(:,6),:);
            Flow(:,ves)=mean(HQFlows)';
            FlowErr(:,ves)=std(HQFlows)';
        end
    end
    if params.PltFlag==1
        plot_MeanVesselFlows(time,Flow,FlowErr)
        if params.SaveData==1
            path2data=string(fullfile(params.data_dir,params.subject));
            if ~exist(path2data, 'dir')
                mkdir(path2data)
            end
            saveas(gcf,fullfile(path2data,'Flow_plot.jpg'))
            close
        end
    end
end