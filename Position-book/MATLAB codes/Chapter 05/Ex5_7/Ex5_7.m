clear all
close all
clc

% Schedule Gain Kalman Filter example
% Based on 
% "Minh Duc Pham, Kay Soon Low, Shu Ting Goh, and Shoushun Chen, 
% “Scheduled Gain EKF for Nano-satellite Attitude Determination System”, 
% IEEE Transactions on Aerospace and Electronic Systems, Vol. 51, Issue 2,
% 22 June 2015, pp 1017-1028, doi: 10.1109/TAES.2014.130204."
% 
% Note:
% Gyro Bias in the paper has been ignore for simplification
% Quarternion vector in this matlab file is q = [q_vector;q_scalar]
%
% Required matlab file:
%    cpm.m, qinverse.m, qmultiply

% Constant Parameter
d2r = pi/180; % deg to radian
wearth = 7.292e-5; % earth rotation rate in rad/sec
wnp = [0;0;wearth]; % nadir pointing rotation
wst = [0;0;0]; % static (or no) rotation 

% Simulation Time Setup
dt = 1; % sampling time as 1 sec
tf = 3600;
t = 0:dt:tf;

% Noise Setup
sign = 1*d2r; 
sigq = 0.1*d2r;

% ------------------------------------------
% True quarternion
% q = [q_vector q_scalar]
qtrue = [-0.3244 -0.6786 -0.6411 0.1522];
% pointing condition
wcond = ones(length(t),1);
wcond(1000:1200) = 0;
wcond(2000:2600) = 0;

for ii = 1:length(t)-1
    if wcond(ii) == 1;
        phitrue = eye(3) -0.5*cpm(wnp)*dt;
        wm(:,ii) = wnp + sign*randn(3,1); % gyro measurement with error
    else
        phitrue = eye(3) -0.5*cpm(wst)*dt;
        wm(:,ii) = wst + sign*randn(3,1); % gyro measurement with error
    end
    
    
    qv = qtrue(ii,1:3)';
    qnext = phitrue*qv;
    qscalar = sqrt(1 - norm(qnext)^2);
    qvec = [qnext' qscalar];
    qvec = qvec/norm(qvec);
    qtrue(ii+1,:) = qvec;
    
    % Output quarternion with noise
    qerr = sigq*randn(3,1);
    qerrv = [qerr;sqrt(1-norm(qerr)^2)];
    qout = qmultiply(qtrue(ii,:)',qerrv);
    qmeas(ii,:) = qout;
end
% ------------------------------------------

% SGEKF models Setup
xe = zeros(1,3); % assume q0 = [0;0;0;1];
phi =  eye(3) -0.5*cpm(wnp)*dt;
C = [eye(3) ];
Q = diag([sign sign sign]).^2;
R = diag([sigq^2 sigq^2 sigq^2]);


% SGEKF Gain Matrix Setup
sqrtden = sqrtm((Q-R*(phi*phi-eye(3)))*(Q-R*(phi*phi-eye(3)))+4*phi*phi*Q*R);
den = 0.5*(sqrtden-(Q-R*(phi*phi-eye(3))))+Q;
num = den + R;
Knp = den*inv(num);

Kst = (sqrtm(Q*Q + 4*Q*R) + Q)*inv(sqrtm(Q*Q + 4*Q*R) + Q + 2*R);

for ii = 1:length(t)-1
   % Generate Estiamted Output
   ye = C*xe(ii,:)';
   qest = [ye;sqrt(1-norm(ye)^2)]; %convert to 4x1 quarternion
   qinv = qinverse(qest);
   dy = qmultiply(qmeas(ii,:)',qinv);
   
   if wcond(ii) == 1
       K = Knp;
   else
       K = Kst;
   end
   
   % Update estimated quarternion/state
   dx = K*dy(1:3);
   dq = [dx;sqrt(1-norm(dx)^2)];
   qxe(ii,:) = qmultiply(dq,qest);
   xe(ii,:) = qxe(ii,1:3);
   
   % Propagation   
   A = eye(3) -0.5*cpm(wm(:,ii))*dt;
   xe(ii+1,:) = A*xe(ii,:)';
   
   % Find the estimated quarternion error for plotting
   qinv = qinverse(qxe(ii,:)');
   qxerr(ii,:) = qmultiply(qtrue(ii,:)',qinv);
end


err = 2*asin(qxerr(:,1:3))/d2r;
qrmse = sum(err.^2,2).^.5/sqrt(3);
plot(t(1:length(t)-1)/60,qrmse)
xlabel('\fontsize{13} Time (min)')
ylabel('\fontsize{13}Quarternion RMSE (deg)')
grid on
ylim([0 1])
xlim([0 60])
set(gca,'fontsize',12)