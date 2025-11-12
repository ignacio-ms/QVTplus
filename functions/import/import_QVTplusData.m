function [data_struct,Labels,v] = import_QVTplusData(path2data)
    DataName = dir(fullfile(path2data,'*.mat'));
    for i=1:length(DataName)
        if length(DataName(i).name)>=7
            if strcmp(DataName(i).name(1:7),'qvtData')
                LOC=i;
                break
            end
        end
    end
    load(fullfile(path2data,DataName(LOC).name));
    clear LOC DataName
    Labels=readLabels(path2data);
    v=Vel_Time_Res;
end