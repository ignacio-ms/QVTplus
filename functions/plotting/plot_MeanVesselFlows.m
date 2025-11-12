function [] = plot_MeanVesselFlows(time,Flow,FlowErr)
    TitleNames={'L ICA';'R ICA';'L MCA';'R MCA';'L ACA';'R ACA';'L PCA';'R PCA';'BA'};
    h1=figure();
    set(h1,'Position',[50 150 750 600])
    order=[1,2,9,3,5,7,4,6,8];
    minY=min(Flow(:));
    maxY=max(Flow(:));
    for ves=1:9
        MEAN=Flow(:,order(ves));
        STD=FlowErr(:,order(ves));
        subplot(3,3,ves)
        h=fill([time';flipud(time')],[MEAN-STD;flipud(MEAN+STD)],'k','Linestyle','None');
        set(h,'facealpha',.3)
        hold on
        plot(time,MEAN,'k')
        xlabel('time (s)','FontSize',8,'FontWeight','Bold')
        ylabel('Flow (mL/s)','FontSize',8,'FontWeight','Bold')
        title(TitleNames{order(ves)},'FontSize',10,'FontWeight','Bold')
        ylim([(minY-0.1.*maxY) (maxY+0.1.*maxY)])
        box off
        grid off
        yticks(0:0.1:1)
        xticks(0:0.15:1.5)
    end
end