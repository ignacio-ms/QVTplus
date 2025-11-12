function [PI_scat,G,XYZ] = compute_length(Mats,BranchList,PI,Area,Quality)
    % compute_length processes the distance of points from their root vessels.
    %
    %       Mats: the oxg connectivity matrix of o outlets, and max (g) generations fromt root to farthers outlet
    % Branchlist: nx6 array with n cross sections, which store [x,y,x,vessel number, vessel cross section location, and global array index [1:n]']
    %         PI: nx1 PI values from n cross sections
    %  PIvel_val: nx1 PI computed from velocity values from n cross sections.
    %             If not computed, leave as PI
    %    Quality: nx1 Quality values from n cross sections.
    %
    % Outputs: [PI_scat,PIvel_scat,G]
    %   PI_scat: cx4 for c points processed from connectivity matrix vessels,
    %       which store [position from root, PI, global cross section number index
    %       (from branch list), and Quality]
    %   PI_scat: same as PI_scat, but storing velocity PI
    %   G: storred a cx1 cell with values as strings of the vessel numer. Can
    %       be used for grouping in plots
    %
    % Used by: PITC_Step1.m
    % Dependencies: None
    
    PI_scat=[]; %Preallocate Empty Array for PI flow
    G = []; %Preallocate Empty Group Array (group is vessel Number)
    for root=1:length(Mats) %for each root
        ConMat=Mats{root}; %Grab the root connectivity
        if ConMat{1,1}>0 %Make sure there is a vessel in connectiviy matrix
            for i=1:length(ConMat(:,1)) %For each row
                for j=1:length(ConMat(1,:)) %For each generation
                    val=ConMat(i,j); %Grab the tuple
                    val=val{:}; %make into a 2D vector
                    if val(2)==2000 %Our vessel flag is always 2000 which means no vessel
        	            val=[];
                    else
                        val=val(2); %assign value 
                    end
                    ConMat(i,j)={val}; %assign the connected vessel
                end
            end
            numGen=0; %number of generations to process foer length
            for j=1:length(ConMat(1,:)) %go through each column
                val=cell2mat(ConMat(:,j)); %get all values in a column
                val(val==2000)=[]; %set 2000 value to empty
                if length(val)>0 %If there are ANY vessels to process in that column, update generation length
                    numGen=j;
                end
            end
            ConMat=ConMat(:,1:numGen); %Trim the ConMat
            for j=1:numGen %for each column
                for i=1:(length(ConMat(:,1))-1) %go through each row 
                    val1num=ConMat{i,j}; %load first value
                    val2num=ConMat{i+1,j}; %load second value
                    if val2num==0 %if we have an empty vessel, it means the vessel is rooted as identical to above
                        ConMat{i+1,j}=val1num; %update that vessel with the above row
                    end
                end
            end
            SegLengths=zeros(size(ConMat)); % Initialise length matrix
            EndP=cell(size(ConMat)); %Make an end point empty cell
            CNT = 1; %count for each new point processed
            for j=1:numGen %For each generation
                for i=1:length(ConMat(:,1)) %for each generation row
                    val1num=ConMat{i,j}; %grab the vessel number
                    if ~isempty(val1num) %If there is a vessel number
                        if SegLengths(i,j)==0 %if the segment value hasn't been processed
                            Data=BranchList(BranchList(:,4)==val1num,:); %Grab all points in vessel
                            %This is only checked after generation 1. It is
                            %used to add extra length connecting the previous
                            %vessels lengths (in prev gen of SegLengths) and
                            %the distance from that vessels end point to the
                            %new vessels start point
                            if j>1 %if Gen 2 or more
                                len=SegLengths(i,j-1); %Initialise length
                                Pend=EndP{i,j-1}; %grab end point of last vessel
                                Pstart=Data(1,1:3); %Grab start point of new vessel
                                Pvec=Pstart-Pend; %get vector of end to start
                                dist=sqrt(Pvec*Pvec'); %compute distance
                                len=len+dist; %add distance to initial seg length
                            else
                                len=0; %This otherwise starts the root with zero length
                            end
                            % Now for each vessel point from Data
                            for k=1:length(Data(:,1))
                                PI_scat(CNT,:,root)=[len , PI(Data(k,6)) , Data(k,6) , Quality(Data(k,6)), Area(Data(k,6))]; %Store length, pulsatility, vessel location index in branchlist, and quality
                                XYZ(CNT,:,root)=[Data(k,1:3)];
                                G{CNT,1,root}=num2str(Data(k,4)); %Store vessel number as a string as well for categorical plotting
                                CNT=CNT+1; %increase count
                                if k < length(Data(:,1))
                                    Pvec=Data(k,1:3)-Data(k+1,1:3); %make point1 to point2 vector
                                    dist=sqrt(Pvec*Pvec'); %Compute distance
                                    len=len+dist; %add distance to length
                                end
                            end
                                EndP{i,j}=Data(k,1:3); %store the end point
                                SegLengths(i,j)=len; %store the length
                        end
                    end
                    if i < length(ConMat(:,1)) %For rows less than end
                        val2num=ConMat{i+1,j}; %Look at the next ConMat value
                        if (val1num-val2num)==0 %If the two vessel roots are the same
                            SegLengths(i+1,j)=SegLengths(i,j); %add the length
                            EndP(i+1,j)=EndP(i,j); %add the end point
                        end
                    end
                end
            end
        else
            PI_scat(1,:,root)=[0 0 0 0 0];
            G{1,1,root}=0;
            XYZ(1,:,root)=[0 0 0];
        end
    end
end