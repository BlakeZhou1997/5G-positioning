clear all
clc
close all
pause(2)

% A simpler version of UAV based Weighted Measurement Fusion KF
% Based on the paper
%          Shu Ting Goh, Ossama Abdelkhalik and Seyed A. (Reza) Zekavat, 
% “A Weighted Measurement Fusion Kalman Filter Implementation for UAV 
% Navigation”, Aerospace Science and Technology, Vol. 28, Issue 1, July 2013,
% pp 315-323, doi: 10.1016/j.ast.2012.11.012
%
% Require matlab files
%     beacone location.m, genmeas, find_delay.m, noise_cov.m, polar2car.m,
%     fusemeas.m

%---------------------------------------------------------------
%Error Parameters
%---------------------------------------------------------------
%Measurement Noise Parameters
d2r=pi/180;
c = 3*10^5;

% Beacon Configuration
range_limit = 30;
beacon_location; % Rbeacon, first row is airport base location


% True Initial Location
rt = [.1 -.1 9.2 .06 .055 -.0048]; %true location when beacon sends out
rta = rt; %true location when UAV receives all beacon signal
%-------------------------------------------
%Total Time for Propagation
%-------------------------------------------
tf=1800;
dt=.1;
t=0:dt:tf;
Tmat = 0; % Time vector for state estimate propagation
m=length(t);

% Noiase Parameter
sigr = 0.1/c; % TOA in seconds
sigp = 0.1*d2r; %DOA
sigv = 1e-3; %process noise
rr = [sigr sigp sigp];

% WMFKF Configuration
inierr = [10*randn(1,3) .01*randn(1,3)]; % Generate random error for initialization
xef = rt(1,:)+inierr; % initial estimate with some random error w.r.t true value
Pf = [15*ones(1,3) .1*ones(1,3)].^2; % initial state erro coveriance
pcovf = diag(Pf);
Qd = sigv^2*[(1/3)*(dt^3)*eye(3) (1/2)*(dt^2)*eye(3); (1/2)*(dt^2)*eye(3) dt*eye(3)]; %process noise

for i=1:m-1
    %-------------------------------------------------------------------
    %Generate measurement
    [ya,tofre,yre]=genmeas(rt(i,:),Rbeacon,c);
    % ya = TOA with signal time delay, and DOAs (elevation and azimuth)
    % tofre = time of flight os ginal
    % yre = beacon range
    
    % We find the beacons that is within given range limit
    ind = find(yre <= range_limit);
    [row,col]=size(ya(:,ind));
    ya_final = reshape(ya(:,ind),row*col,1); %doing some reshaping
    
    %we only include the beacon within 10km radius
    err =kron(rr,ones(length(ind),1)).*randn(length(ind),3);
    err = reshape(err',row*col,1);
    ym = ya_final+err;
    tof = tofre(ind); %use for measurement fusion weight calculator
    %only beacon position information has noise
    Rbwnoise = Rbeacon(2:end,:) + .001*randn(24,3);
    Rbmeas = [Rbeacon(1,:);Rbwnoise];
    Rbmeas = Rbmeas(ind,:);    
    %-------------------------------------------------------------------
    
    % ---------------------------------------
    % True Model Propagation        
    A = [eye(3) dt*eye(3);zeros(3) eye(3)];
    rt(i+1,:)= (A*rt(i,:)')';    
    
    dtdelay = dt + max(tof);    
    A = [eye(3) dtdelay*eye(3);zeros(3) eye(3)];
    rta(i+1,:)= (A*rt(i,:)')';
    % ---------------------------------------    
    
    % -----------------------------------------------------
    %Fusion method
    %Based on Last DOA and in cartesian
    % ------------------------------------------
    %find the furthest beacon
    [~,maxtof(i,:)]= max(tof);
    Tmat(i+1,:) = t(i) + max(tof);
    dtdelay = Tmat(i+1,:) - Tmat(i,:);
    A = [eye(3) dtdelay*eye(3);zeros(3) eye(3)];
    
    %Propagate to the time where last measurement received
    xef(i+1,:) = (A*xef(i,:)')';
    pcovf = A*pcovf*A'+Qd;
    Pf(i+1,:)=diag(pcovf);
    
    % Do measurement fusion
    [yf,Rf,H,~,~]=fusemeas(ym,tof,Rbmeas,sigr*c/2,sigp,c,maxtof(i,:));
    yef = Rbmeas(maxtof(i,:),:)'-xef(i+1,1:3)';
    
    gain = pcovf*H'*inv(H*pcovf*H'+Rf);
    pcovf = (eye(6)-gain*H)*pcovf;
    xef(i+1,:)=xef(i+1,:)+(gain*(yf-yef))';
    
end



% Three sigma axis
sig3 = 1000*3*Pf.^.5;

% Plot error w.r.t to true location when UAV receives all beacon signal
figure(1)
subplot(311)
plot(t/60,1000*xef(:,1)-1000*rta(:,1),t/60,-sig3(:,1),t/60,sig3(:,1))
xlabel('\fontsize{13}Time (min)')
ylabel('\fontsize{13}x-err (m)')
set(gca,'FontSize',12)
axis([0 30 -60 60])
subplot(312)
plot(t/60,1000*xef(:,2)-1000*rta(:,2),t/60,-sig3(:,2),t/60,sig3(:,2))
xlabel('\fontsize{13}Time (min)')
ylabel('\fontsize{13}y-err (m)')
set(gca,'FontSize',12)
axis([0 30 -60 60])
subplot(313)
plot(t/60,1000*xef(:,3)-1000*rta(:,3),t/60,-sig3(:,3),t/60,sig3(:,3))
xlabel('\fontsize{13}Time (min)')
ylabel('\fontsize{13}z-err (m)')
set(gca,'FontSize',12)
axis([0 30 -60 60])

