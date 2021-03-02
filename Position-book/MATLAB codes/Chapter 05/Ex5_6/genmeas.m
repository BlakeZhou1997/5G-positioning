function [ym,TOA,yre]=genmeas(rt,Rm,c)

% Generate all TOA and DOA measurements w.r.t to all available beacons
% measurements include the TOF of signal and the moment UAV receives the
% signal
% Output
%    ym - containt TOA and DOA
%    TOA - time of flight of signals
%    yre - beacon range to UAV
%
%
% Required matlab file
%      find_delay.m

[row,~]=size(Rm);
TOF1 = 1/c*sum((kron(rt(1:3),ones(row,1))-Rm).^2,2).^.5;


for jj=1:row
    dtp = (TOF1(jj));
    A = [eye(3) dtp*eye(3);zeros(3,6)];
    rt3 = (A*rt')';    
    TOF2(jj,1) = find_delay([Rm(jj,1:3)-rt3(1:3) rt3(4:6)],c);
    dtp = TOF2(jj,1);
    A = [eye(3) dtp*eye(3);zeros(3,6)];
    rt4 = (A*rt3')';

    DOA1(jj,1) = atan2(Rm(jj,2)-rt4(2),Rm(jj,1)-rt4(1));
    DOA2(jj,1) = atan2(Rm(jj,3)-rt4(3),sqrt((Rm(jj,1)-rt4(1))^2+(Rm(jj,2)-rt4(2))^2));
    yre(jj,:) = TOF2(jj,1)*c;   
   
end


 

TOA = TOF1+TOF2;
DOA1 = DOA1;
DOA2 = DOA2;
ym = [TOA';DOA1';DOA2'];

