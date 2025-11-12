function Mats = compute_conn(startroots,BranchList,SearchDist,LR)
% compute_conn creates the connectivity matrix for floating vessels
% processed through QVT.
%
% startroots: 1xn array of vessel numbers to computer connectivity from 
% BranchList: nx6 array with n cross sections, which store [x,y,x,vessel number, vessel cross section location, and global array index [1:n]']
% SearchDist: Search distance for connectivity
%         LR: array of left and right ordered vessel numbers. Leave [] if empty.
%
% Outputs: [Mats]
%   Mats: 1xn cell, where n is number of root vessels, and each cell
%   contains a oxg matrix of o outlets, and g  vessel generations to get to all
%   outlets from root vessel
%
% Used by: PITC_Step1.m
% Dependencies: search_local_points.m, pick_point_crossover.m, insertrows.m
    %% load start and end vessels and compute tangents
    blist=unique(BranchList(:,4)); %Branch List with all info positions, vessel seg number, 
    numbranches=length(blist); %number of branches
    Terms=zeros([2*numbranches 8]); %Number of Start and End Points
    for i=1:numbranches %For each branch
        Data = [BranchList(BranchList(:,4)==blist(i),:)]; %Collect all points
        P1s = Data(1,1:3); %Start point
        P1e = Data(3,1:3); %two points up from start
        T1 = P1e-P1s; %Start Point Tangent
        t1 = T1./(sqrt(T1*T1')); %start point unit tangent
        P2s = Data((end-2),1:3); %end point start (2 points back)
        P2e = Data(end,1:3); %end point end
        T2 = P2e-P2s; %End point tangent
        n2 = T2./(sqrt(T2*T2')); %end point unit tangent
        pos=1+(i-1)*2; %Index Location in Terms (goes up by two)
        Terms(pos,:) = [blist(i),0,P1s,t1]; %Terms ordered vessel num, 0/1=start/end, Point location, tangent
        Terms(pos+1,:) = [blist(i),1,P2e,n2];
    end
    %% Compute Connectivity for each root
    Avoid=[]; %Avoid vector stores any vessels that are processed so we don't double dip
    flg=0;
    for root=1:length(startroots) %ICA's and basilar
        startr = startroots{root}; %Start at whatever root based on startroot vessel numbers
        ConVec = cell([1,20]); %Empty connectivity matrix with up to 20 generations (no more than than will ever be needed)
        ConVec(:,:) = {[0 0]}; %Initialise all cells (this will help identify filled slots later)
        numIter=0; %Number of iterations (generations)
        if ~isempty(startr) %As long as there IS a startroot vessel
            for j=1:length(startr) %In case startroot is a collection of vessels
                ConVec(1,j)={[startr(j),startr(j)]}; %Make a connection
                start=startr(j); %Whatever final vessel this stores will start the real connectivity code next
                numIter=1+numIter; % Increase generation (column)
            end
            Flag=1; %Keep the loop running until there are no new connectivities
            while Flag == 1
                cnt=0; %Count of how many new vessels are added this generation
                numIter=numIter+1; %Start new generation
                newstart=[]; %Initialise the vector of new start roots for the next generation
                for i=1:length(start) %go through start vectors
                    Pend=Terms(Terms(:,1)==start(i) & Terms(:,2)==1,3:5); %End point of vessel (where we will look for new ones)
                    %P_options=[]; %Initialise Point options matrix
                    Avoidnew=unique(cell2mat(ConVec)); %Get new vessels to avoid based on what was connected last generation
                    [aa,bb]=size(Avoidnew);
                    if aa<bb
                        Avoidnew=Avoidnew'; %Make sure we have a column matrix
                    end
                    Avoid=[Avoid;Avoidnew]; %Keep just avoided vessels to date
                    Avoid=unique(Avoid); %Keep avoid vector slim
                    if sum(size(LR))>0 %add on Left and Right avoids at last step
                        if root==1 %Left ICA
                            AvoidF=[Avoid;LR(:,2)]; %avoid right vessels
                        else
                            AvoidF=[Avoid;LR(:,1)]; %avoid left vessels
                        end
                    else
                        AvoidF=Avoid; %If there is no Left Right Difference
                    end
                    Options=Terms; %Load Options as all potential points
                    for j=1:length(AvoidF)
                        [preused,~]=find(Options(:,1)==AvoidF(j)); %Find index of points to avoid
                        if ~isempty(preused)
                            Options(preused,:)=[]; %Empty out the options
                        end
                    end
                    if flg==0
                        if root==3
                            SearchDist=SearchDist*1.0;
                            flg=1;
                        end
                    end
                    [Options]=search_local_points(Pend,Options,SearchDist); %Return point options from Options around Pend within SearchDist
                    [r,~]=find(Options(:,2)==1); %Look for options that are terminals
                    if length(r)>1 % More than 1 terminal, means a crossover
                        [r,~]=find(Options(:,2)==0); %Find the start points
                        OptionsN=Options(r,:); %OptionsN are the real vessel options
                        if ~isempty(r) %If there aren't any start points (0), then it's two terminals pointing at each other
                            Pstart=Terms(Terms(:,1)==start(i) & Terms(:,2)==1,:); %load the entire end point with tangent
                            P_options = pick_point_crossover(Pstart,OptionsN); %from options, Pick the ideal crossover
                        else
                            P_options = Options([],:); %No actual start vessel options
                        end
                    else
                        [r,~] = find(Options(:,2)==0); %Grab all start vessels- indexs (avoining initial terminal)
                        P_options = Options(r,:); %All start vessels are the options
                    end
                    % End of the day, We have P_options and the rest of this can be followed smoothly
                    numBif=length(P_options(:,1));
                    Prev=cell2mat(ConVec(:,numIter-1)); % look at previous column of connections
                    Prev=Prev(:,2); %get the right half of connection tuple
                    [idx]=find(Prev==start(i)); %Find the row where index where the connection is to insert below
                    if numBif>1 %insert rows below for new connections off this branch
                        [ConVec]=insertrows(ConVec,numBif,idx);
                    end  
                    for Bifs=1:numBif
                        ConVec(idx-1+Bifs,numIter)={[start(i) P_options(Bifs,1)]};
                        cnt=cnt+1;
                    end
                    if ~isempty(P_options)==1
                        newstart=[newstart;P_options(:,1)]; %Add on newly connected vessels to the new start vector
                    else
                        ConVec(idx,numIter:end)={[-1,2000]}; % CHANGED FROM 1 If a vessel was not connected to anythign new, fill rest of connectivity row with filler value
                    end
                end
                start=newstart; %Update new  start
                if cnt>0     
                    Flag=1; %If new vessels have been addded, keep it going
                else
                    Flag=0; %Else, stop the loop
                end
            end
        end
        Mats(root)={ConVec}; %Save the root connectivity vector and start again
    end
end