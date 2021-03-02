clear all
close all
clc

% Example 19.2
% Case 2 Scenario 1 from the following paper
% Shu Ting Goh, Ossama Abdelkhalik and Seyed A. (Reza) Zekavat, “Spacecraft 
% Constellation Orbit Estimation via a Novel Wireless Positioning System”, 
% 2009 AAS/AIAA Space Flight Mechanics Meeting, Feb 2009, pp 1-19, 
% Savannah, Georgia, AAS 09-116
%
%
% Require matlab files:
%     basis_function.m, orbit.m, covprop.m, diff_orbit.m, doa2.m, orb2eci.m

format long
mu=398600.4; % earth gravitational constant
d2r = pi/180; % degree to radian

%-------------------------------------------
%Total Time for Propagation
%-------------------------------------------
tf = 12*3600; % 12 hours of simulation
dt = 10; % sampling rate at 10 seconds
t=0:dt:tf; % create a time vector
m=length(t);

%---------------------------------------------------------------
%Error Parameters
%---------------------------------------------------------------
% Measurement Noise Parameters
% WLPS noise
sigtoa = 0.01; % 10 meter
sigdoa = 0.01*pi/180; % 0.01 degree
% Ground Radar noise
siggtoa = 0.001; % 1 meter
siggdoa = .001*pi/180;
yerr = kron([siggtoa siggdoa siggdoa sigtoa sigdoa sigdoa],ones(m,1)).*randn(m,6); % error matrix for ground radar

% Process Noise Parameters
qx = 1e-6;
qy = qx;
qz = qx;


%-------------------------------------------------------------------
%Model Initial Condition
%-------------------------------------------------------------------
% Setup spacecraft 1 and 2
oev1 = [12000+6378.145 0.3 25*d2r 30*d2r 20*d2r 30*d2r];
[r1,v1] = orb2eci(oev1);
oev2 = [12000+6378.145 0.3 25*d2r 30.1*d2r 20*d2r 30*d2r];
[r2,v2] = orb2eci(oev2);
sc1 = [r1;v1]';
sc2 = [r2;v2]';

% Kalman filter initialization
% initial state vector and error covariance
xe = [r1;v1;r2;v2]';
pcov = [1 1 1 0.5 0.5 0.5 1 1 1 0.5 0.5 0.5].^2;
% covariance for measurement and process noise
P = diag(pcov);
R = diag([siggtoa siggdoa siggdoa sigtoa sigdoa sigdoa]).^2;
Q = diag([qx qy qz qx qy qz].^2);
D = [zeros(3,6);
    eye(3) zeros(3);    
    zeros(3,6);
    zeros(3) eye(3);];


for i = 1:m-1
    
    %----------------------------------------------------------------
    %Ground radar for SC1, TOA, DOA between sc1 and sc2
    %----------------------------------------------------------------
    % --- Generate True Measurement ---
    % Ground Radar
    r1 = norm(sc1(i,1:3));
    doa1 = doa2(sc1(i,1:3)');
    
    % WLPS
    r12 = sc2(i,1:3) - sc1(i,1:3);    
    toa12 = norm(r12);
    doa12 = doa2(r12);    
    
    % Measurement vector
    ym = [r1;doa1;toa12;doa12] + yerr(i,:)';
    
    %----------------------------------------------------------------
    % Kalman filter
    %----------------------------------------------------------------
    % Generate estimated output
    % Est. Ground Radar
    r1e = norm(xe(i,1:3));
    doa1e = doa2(xe(i,1:3)');    
    % WLPS
    r12e = xe(i,7:9) - xe(i,1:3);    
    toa12e = norm(r12e);
    doa12e = doa2(r12e);      
    %Estimate Output
    ye = [r1e;doa1e;toa12e;doa12e];
    
    
    %Basis function/Sensitivity Matrix
    C=basis_function(xe(i,:)');   
    
    %Gain, Covariance and States update
    K = P*C'*inv(C*P*C'+R);
    P = (eye(12) - K*C)*P;
    xe(i,:) = xe(i,:)+(K*(ym-ye))';
       
    %------------------------------------------
    %Estimated State Propagation
    %------------------------------------------
    xe(i+1,:) = zeros(1,12);
    x0 = xe(i,1:6)';
    [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
    xe(i+1,1:6) = xx(end,:);
    x0 = xe(i,7:12)';
    [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
    xe(i+1,7:12) = xx(end,:);
    
    %------------------------------------------
    %Covariance Propagation in Time-Continuous
    %------------------------------------------    
    F = [zeros(3) eye(3) zeros(3,6);
        diff_orbit(xe(i,1:3),mu) zeros(3,9);
        zeros(3,9) eye(3);
        zeros(3,6) diff_orbit(xe(i,7:9),mu) zeros(3)];
    
    p0 = reshape(P,12*12,1);
    [~,xx] = ode45(@(tt,pl)covprop(tt,pl,F,D,Q),[0 dt],p0);
    pfinal = xx(end,:);
    P = reshape(pfinal,12,12);
    pcov(i+1,:) = diag(P);   
    
    %-------------------------------------------------------------------
    %Truth States Propagation 
    %-------------------------------------------------------------------
    x0 = sc1(i,:)';
    [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
    sc1(i+1,:) = xx(end,:);
    x0 = sc2(i,:)';
    [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
    sc2(i+1,:) = xx(end,:);
end

% Calculate estimation error and do plotting
sc1err = sc1 - xe(:,1:6);
sc2err = sc2 - xe(:,7:12);

p1rmse = sum(sc1err(:,1:3).^2,2).^.5/sqrt(3);
p2rmse = sum(sc2err(:,1:3).^2,2).^.5/sqrt(3);

figure(1)
plot(t/60,p1rmse*1000,t/60,p2rmse*1000,'-.')
xlabel('\fontsize{13}Time (min)')
ylabel('\fontsize{13}Position RMSE (m)')
set(gca,'fontsize',12)
legend('S/C 1','S/C 2')
grid on
ylim([0 100])
xlim([0 720])

figure(2)
subplot(311)
plot(t/60,sc1err(:,1)*1000,t/60,3*pcov(:,1).^.5*1000,t/60,-3*pcov(:,1).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-x error (m)')
ylim([-100 100])
xlim([0 720])
grid on
subplot(312)
plot(t/60,sc1err(:,2)*1000,t/60,3*pcov(:,2).^.5*1000,t/60,-3*pcov(:,2).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-y error (m)')
ylim([-100 100])
xlim([0 720])
grid on
subplot(313)
plot(t/60,sc1err(:,3)*1000,t/60,3*pcov(:,3).^.5*1000,t/60,-3*pcov(:,3).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-z error (m)')
ylim([-100 100])
xlim([0 720])
grid on

figure(3)
subplot(311)
plot(t/60,sc2err(:,1)*1000,t/60,3*pcov(:,7).^.5*1000,t/60,-3*pcov(:,7).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-x error (m)')
ylim([-100 100])
xlim([0 720])
grid on
subplot(312)
plot(t/60,sc2err(:,2)*1000,t/60,3*pcov(:,8).^.5*1000,t/60,-3*pcov(:,8).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-y error (m)')
ylim([-100 100])
xlim([0 720])
grid on
subplot(313)
plot(t/60,sc2err(:,3)*1000,t/60,3*pcov(:,9).^.5*1000,t/60,-3*pcov(:,9).^.5*1000)
xlabel('\fontsize{12}Time (min)')
ylabel('\fontsize{12}pos-z error (m)')
ylim([-100 100])
xlim([0 720])
grid on