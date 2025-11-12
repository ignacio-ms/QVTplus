function [P_options]=search_local_points(Pend,Terms,SearchDist)
% search_local_points finds any floating points that may continue the
% vessel from Pend
%
% Pend: 1x6 vector with data [x,y,z,x',y',z']
% Terms: nx8 matrix if n potential points with data [global vessel location, flag for inlet or outlet (1/0), and x,y,z, and x',y',z']
% SearchDist: Voxel search box distances for points around Pend
%
% Outputs: [P_options]
%   P_options: the single new vessel that continues Pstart
%
% Used by: compute_conn.m
% Dependencies: None
    [r,~]=find(Terms(:,3)<(Pend(1)+SearchDist));
    TempPs=Terms(r,:);
    [r,~]=find(TempPs(:,3)>(Pend(1)-SearchDist));
    TempPs=TempPs(r,:);
    [r,~]=find(TempPs(:,4)<(Pend(2)+SearchDist));
    TempPs=TempPs(r,:);
    [r,~]=find(TempPs(:,4)>(Pend(2)-SearchDist));
    TempPs=TempPs(r,:);
    [r,~]=find(TempPs(:,5)<(Pend(3)+SearchDist));
    TempPs=TempPs(r,:);
    [r,~]=find(TempPs(:,5)>(Pend(3)-SearchDist));
    P_options=TempPs(r,:);
end