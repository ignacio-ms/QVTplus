function RealVal = readLabels(path)
% readLabels reads the saved input table from QVT+,
%
% path: All that's needed is the path to the "LabelsQVT.csv"
%
% Outputs: [RealVal]
%   RealVal: A poorly named matrix which is just the cells of the original
%   table after being read from the csvs
%
% Used by: autoCollectFlow.m, and any separate functions to compute any
% saved data (PITC codes, Damping codes etc)
%
% Dependencies: None    
    if exist(fullfile(path,'LabelsQVT.csv')) == 2
        fprintf('Loaded previous labels.\n')
        cont=1;
    else
        fprintf('Labels have not been processed for this subject.\n')
        cont=0;
    end
    %% This function just parses the filled out labelledbranches form for all branch numbers
    %% It reads line bny line then outputs the branch numbers and names into a cell
    if cont==1
        fid = fopen(fullfile(path,'LabelsQVT.csv'));
        tline = fgetl(fid);
        RealVal=cell([10,3]);
        for i=1:10
            tline = fgetl(fid);
            A = strsplit(tline,',');
            LOC=1;
            MergeArray=[];
            for j=1:length(A)
                mark1=strsplit(A{j},'[');
                mark2=strsplit(A{j},']');
                if length(mark1)>1
                    MergeArray(LOC,1)=j;
    
                end
                if length(mark2)>1
                    MergeArray(LOC,2)=j;
                    LOC=LOC+1;
                end
            end
            B=cell([1 3]);
            B(1)=A(1);
            B(2)={erase(A{2},"""")};
            B(3)={erase(A{end},"""")};
            if length(MergeArray)>0
                MergeLen=length(MergeArray(:,1));
                for j=1:MergeLen
                    if MergeArray(j,1)==2
                        value='';
                        for k=MergeArray(j,1):(MergeArray(j,2)-1)
                            value=strcat(value,A{k},',');
                        end
                        value=strcat(value,A{k+1});
                        value = erase(value,"""");
                        B(2)={value};
                    else
                        value='';
                        for k=MergeArray(j,1):(MergeArray(j,2)-1)
                            value=strcat(value,A{k},',');
                        end
                        value=strcat(value,A{k+1});
                        value = erase(value,"""");
                        B(3)={value};
                    end
                end
            end
            RealVal(i,:)=B;
        end
        fclose(fid);
    end
end