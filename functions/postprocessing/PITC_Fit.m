% PITC_Fit fits the pulsatiltiy length PITC slope. It uses quality cutoff. 
%
% All that's needed is the path to QVT saved data which must include
% the RawDamping.mat filename (should be renamed(
%
% Outputs: The PITC slopes for each root, it can optionally save a figure
% with all 3 fits.
plotFlag=1; %keep 1 if you want to plot and save results to population
path2data='C:\Users\sdem348\Desktop\Gonzalo\derivatives\QVT\sub-001\breathing'; % to the QVT data folder
%% Crop Data by Quality>Thresh as a Baseline (Plot as Well)
load(fullfile(path2data,'RawDamping.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Actual Processing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
ResCell=cell([1 3]); %This is what we store thresholded values into
thresh=2.5; %Quality Threshold
for i=1:3 %three roots
    [idx,~]=find(PI_scat(:,4,i)>thresh);
    Vals=squeeze(PI_scat(idx,:,i));
    ResCell(i)={abs(Vals)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if plotFlag==1 %Optional
    MrkrSz=7; %Plotting MarkerSize
    h1=figure(1);
    set(h1,'Position',[50 150 1200 600])
    TitleNames={'Left ICA';'Right ICA';'Basilar'};
    for i=1:3 %three roots
        Vals=squeeze(PI_scat(:,:,i));
        subplot(2,3,i)
        scatter(Vals(:,1),Vals(:,4),MrkrSz,'MarkerFaceColor','k','MarkerEdgeColor','k',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0)
        xticks(0:25:400)
        yticks(0:.5:4)
        xlim([0 max(Vals(:,1))])
        ylim([0 4])
        Vals=ResCell{i}; %Now Use Thresholded Values
        hold on
        scatter(Vals(:,1),Vals(:,4),MrkrSz,'MarkerFaceColor','b',...
        'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
        xlabel('distance from root (mm)','FontSize',10,'FontWeight','Bold')
        title(TitleNames{i},'FontSize',12,'FontWeight','Bold')
        if i==1
            ylabel('Quality (max 4)','FontSize',10,'FontWeight','Bold')
            L1=legend('Raw',strcat('Thresholded > ',num2str(thresh)),'Location','SouthWest');
            set(L1,'FontSize',8)
        end
        grid on
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if plotFlag==1 %Optional
    MrkrSz=7; %Plotting MarkerSize
    TitleNames={'Left ICA';'Right ICA';'Basilar'};
    for i=1:3 %three roots
        Vals=ResCell{i};
        subplot(2,3,i+3)
        scatter(Vals(:,1),Vals(:,2),MrkrSz,'MarkerFaceColor','b',...
        'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
        xticks(0:25:400)
        yticks(0:.2:2)
        xLims=[0 max(Vals(:,1))];
        xlim(xLims)
        ylim([0 1.3])
        %if max(Vals(:,2))>1
        %    ylim([0 1.8])% 1.05*max(Vals(:,2))])
        %else
        %    ylim([0 1])
        %end
        hold on
        xlabel('distance from root (mm)','FontSize',10,'FontWeight','Bold')
        title(TitleNames{i},'FontSize',12,'FontWeight','Bold')
        if i==1
            ylabel('PI_{flow}','FontSize',10,'FontWeight','Bold')
            lgnd=legend(strcat('Quality >',num2str(thresh),' data'),'');
            set(lgnd,'Location','SouthWest')
        end
        grid on
    end
end
%% Fit Slopes (Linear)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Actual Processing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
ResSlope=cell([1 3]);
STORE=[];
for i=1:3
    Vals=ResCell{i};
    posx=Vals(:,1);
    posy=Vals(:,2);
    W=(Vals(:,4)-thresh)./(4-thresh);
    meanPI=sum(W.*Vals(:,2))./sum(W);
    STORE=[STORE;W Vals(:,2)];
    PIs(1,i)=sum(W.*Vals(:,2))./sum(W);
    f=fit(posx,posy,'poly1','Weights',W);
    ci = confint(f,0.68);
    C = coeffvalues(f); % Get the coefficient values
    yCalc=f(posx);
    Rsq1 = 1 - sum((posy - yCalc).^2)/sum((posy - mean(posy)).^2);
    ResSlope(i)={[C(1) C(2) Rsq1 ci(1,1) ci(2,1) meanPI]};
end
MPI=sum(STORE(:,1).*STORE(:,2))./sum(STORE(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if plotFlag==1 %Optional
    for i=1:3 %three roots
        Vals=ResCell{i};
        xLims=[0 max(Vals(:,1))];
        Vals=ResSlope{i};
        ypos=Vals(1).*xLims+Vals(2);
        subplot(2,3,i+3)
        hold on
        plot(xLims,ypos,'k-')
        lgnd.String(2)={'PI(dist)=PI_{TC}dist+b'};
        lgnd.AutoUpdate='off';
        set(lgnd,'FontSize',8,'Location','NorthWest')
    end
    saveas(gcf,fullfile(path2data,'PITC_Plot.png'))
end
PITC={ResSlope(1) ResSlope(2) ResSlope(3)};
save(fullfile(path2data,'PITC'),'PITC')



