% Just contracts to values that have measurements and stores location
function [X1n,X2n,loc] = enc_matchSamples(X1,X2)
    count=1;
    X1n=[];
    X2n=[];
    loc=[];
    for i=1:length(X1)
        if X1(i)~=0 && X2(i)~=0
            if ~isnan(X1(i)) && ~isnan(X2(i))
                X1n(count,1)=X1(i);
                X2n(count,1)=X2(i);
                loc(count,1)=i;
                count=count+1;
            end
        end
    end
end