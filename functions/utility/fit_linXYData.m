function [mdl,plotData] = fit_linXYData(xData,yData,W)
if isempty(W)
        mdl = fitlm(xData,yData,'linear');
        rng=[min(xData):(max(xData)-min(xData))/50:max(xData)]';
        [Ymain,yci] = predict(mdl,rng);
        plotData=[rng,Ymain,yci];
else
        mdl = fitlm(xData,yData,'linear','Weights',W);
        rng=[min(xData):(max(xData)-min(xData))/50:max(xData)]';
        [Ymain,yci] = predict(mdl,rng);
        plotData=[rng,Ymain,yci];
end