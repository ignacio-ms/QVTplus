function P_options = pick_point_crossover(Pstart,OptionsN)
% pick_point_crossover picks the most likely point in the point options
% which continues Pstart at crossovers
%
% Pstart: 1x5 vector with data [global vessel location, flag for inlet or outlet (1/0), and x,y,z]
% OptionsN: nx8 matrix if n potential points with data [global vessel location, flag for inlet or outlet (1/0), and x,y,z, and x',y',z']

%
% Outputs: [P_options]
%   P_options: the single new vessel that continues Pstart
%
% Used by: compute_conn.m
% Dependencies: None
    [ps,~]=size(OptionsN);
    N1=Pstart(6:8);
    Likelihood=[];
    for j=1:ps
        Np = OptionsN(j,6:8);
        Pb = OptionsN(j,3:5)-Pstart(3:5);%vector bridging the terminal to condidate point
        Pbn = Pb./(sqrt(Pb*Pb'));%normalized pointing vector 
        %so we have 3 tangents from selected point, the terminal
        %tangent, the condidate point tanged (from continuing
        %along) and the bridge tangent. We will dot the terminal
        %with condidate, and bridge, then take the mean for each
        %candidate, the max mean is the likely continuation of the
        %vessel.
        dot1=N1*Np';
        dot2=N1*Pbn';
        Likelihood(j,:)=[j mean([dot1 dot2])];
    end
    [~,b]=max(Likelihood(:,2));
    P_options=OptionsN(b,:);
end
