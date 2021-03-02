function dbs_pos = gen_dbs_pos(N,Ny,L)
%N = number of antenna area in each circle
N = N+1;
%assume 10 in each each, so 11 for linspace as 360 deg = 0 deg
R = 60; %distance from center of body, in cm


z = linspace(-60,L,Ny);

% ini_theta = linspace(0,0,Ny+1);
ini_theta = linspace(0,90,Ny+1);
% ini_theta = linspace(0,180,Ny+1);

m = length(z);

dbs_pos = [];

for i=1:m
    theta = linspace(0+ini_theta(i),360+ini_theta(i),N);
    theta = theta(1:end-1)*pi/180;
    for kk = 1:N-1
    pos(kk,:) = [R*cos(theta(kk)) R*sin(theta(kk)) z(i)];
    end
    dbs_pos = [dbs_pos;pos];
end


