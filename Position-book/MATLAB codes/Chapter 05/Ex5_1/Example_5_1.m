clear all
clc
close all

% Example 19.1
% This is a matlab file for simple winding rotation tracking


d2r = pi/180; % deg to radian
dt = 1; % Sampling time
xv = 0:30000;

% --- True State ----
xtrue = 5*d2r; % true state location at initial
wtrue = 19*360/60*d2r; %actual rotation speed, 19 RPM

% --- Initialization ----
x = [0 0];
west = 18*360/60*d2r; % estimated/expected rotation speed input, 18 RPM 
% initial state error covariance
p = [10 5]*d2r;
P = diag(p.^2);

% --- State model -----
% dynamic model
A = [1 dt;0 1];
B = [1;0]; % input matrix
% Process noise model
wnoise = 0.01*d2r; 
Q = wnoise^2; % process noise covariance
D = [1;0];
% Measurement model
C = [1 0];
v = 1*d2r; % measurement noise
R = v^2; % measurement nosie covariance

% Estimation up to 1000 steps
for ii = 1:length(xv)-1
    % Generate measurements from the true value
    xtrue(ii+1,:) = mod(xtrue(ii,:) + wtrue*dt,2*pi);
    % mod is used to ensure always within 0 to 360 degrees
    ym = xtrue(ii+1,:) + v*randn(1);
    
    % Estimation model for Kalman filter process
    x(ii+1,:) = (A*x(ii,:)' + B*west*dt)';
    x(ii+1,1) = mod(x(ii+1,1),2*pi);
    % mod is used to ensure always within 0 to 360 degrees
    % Estimated output
    ye = C*x(ii+1,:)';
    
    % Calculate Kalman gain
    K = P*C'*inv(C*P*C'+R);
    
    % Kalman filter update for state vector
    x(ii+1,:) = x(ii+1,:) + (K*(ym-ye))';    
    % Covariance matrix update
    P = (eye(2) - K*C)*P;
    
    % Covariance matrix propagation
    P = A*P*A'+D*Q*D';
    % save diagonal elements of covariance at each step
    p(ii+1,:) = diag(P.^.5);
end

% plot error
figure(1)
plot(xv/60,(xtrue(:,1) - x(:,1))/d2r,xv/60,3*p(:,1)/d2r,xv/60,-3*p(:,1)/d2r)
ylim([-1 1])
grid on
xlabel('\fontsize{13}Time steps (min)')
ylabel('\fontsize{13}Winding location error (deg)')
set(gca,'fontsize',12)

figure(2)
plot(xv/60,(wtrue - west - x(:,2))/d2r,xv/60,3*p(:,2)/d2r,xv/60,-3*p(:,2)/d2r)
ylim([-0.01 0.01])
grid on
xlabel('\fontsize{13}Time steps (min)')
ylabel('\fontsize{13}Rotation bias estimation error (deg/sec)')
set(gca,'fontsize',12)

figure(3)
plot(xv/60,(wtrue - west)/d2r*ones(1,length(xv/60)),'r-.',xv/60,x(:,2)/d2r)
xlabel('\fontsize{13}Time steps (min)')
ylabel('\fontsize{13}Rotation bias estimation (deg/sec)')
grid on
legend('True Bias', 'Estimated Bias')
set(gca,'fontsize',12)