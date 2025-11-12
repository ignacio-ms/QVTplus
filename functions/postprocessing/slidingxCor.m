function [rho]=slidingxCor(csaFlowINT,defCurve)
    newCurve=csaFlowINT;
    rho=zeros([(length(newCurve)-1) 1]);
    for i=1:(length(csaFlowINT)-1)
        newCurve=[newCurve(2:end) newCurve(1)];
        rho(i,1) = corr(newCurve',defCurve');
    end

end
