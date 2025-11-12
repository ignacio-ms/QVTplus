function [ConVec]=insertrows(ConVec,numBif,idx)
% insertrows adds new rows to the connectivity matrix when vessels
% bifurcate, letting connectivity matrix not need to be preallocated any
% size
%
% ConVec: oxg matrix of 0 outlets, and g max vessel generations
% numBif: number of bifurcations at the vessel crossing
%    idx: initial row (outlet) that will be bifucation i.e. position to add
%         new rows
%
% Outputs: [ConVec]
%   ConVec: (o+numBif-1)xg cell, where new blank rows are initialised into the connectivity matrix
%
% Used by: compute_conn.m
% Dependencies: None
    [Rows,~]=size(ConVec);
    Test=cell([(Rows+numBif-1),20]);
    M1=ConVec(1:idx,:);
    M2=ConVec(idx+1:end,:);
    temp=cell([(numBif-1) 20]);
    temp(:,:)={[0 0]};
    Test(:,:)=[M1;temp;M2];
    ConVec=Test;
end