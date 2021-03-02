clear all
close all
clc

% Constant Parameter
d2r = pi/180; % deg to radian
f = 4.0*10^9; %c-band communication frequency
c = 3*10^5;

% Simulation time configuration
dt = 1;
tf = 300;
t = 0:dt:tf;

% Vehicle Parameter and Dynamic Model
towerloc = [20 25;
             -5 25]; % communication tower location
xtrue = [30 50 -0.1 -0.12]; % vehicle initial location
u = [-0.25e-3;0]; % acceleration
B = [0.5*dt^2*eye(2);dt*eye(2);];
A = [eye(2) dt*eye(2);
    zeros(2) eye(2);];

% Noise parameter
sigd = 0.5*d2r; % doa error
sigf = 10; % doppler shift error
sigv = [sigf;sigd;sigd];
siga = 1e-5; % accelerometer noise, 10cm/s^2
R = diag(sigv).^2;

% Ensemble Kalman filter
N = 64;
inierr = [3;3;0.1;0.1]; % Initial error or distribution
X = kron(xtrue',ones(1,N)) + kron(inierr,ones(1,N)).*randn(4,N);

% EKF
iniekf = [1;1;0.1;0.1];
xekf = xtrue + iniekf'.*randn(1,4);
Pekf = diag(iniekf.^2);
pekf = diag(Pekf);
Qd = 1e-4*[(1/3)*(dt^3)*eye(2) (1/2)*(dt^2)*eye(2); (1/2)*(dt^2)*eye(2) dt*eye(2)]; %process noise


for ii = 1:length(t)-1
    
    % Generate measurement    
    dopfreq = norm(xtrue(ii,3:4))/c*f;
    r1 = xtrue(ii,1:2) - towerloc(1,:);
    r2 = xtrue(ii,1:2) - towerloc(2,:);
    ym = [dopfreq;doa(r1);doa(r2)] + sigv.*randn(length(sigv),1);    
    
    %% Ensemble Kalman filter
    % Add oscillation into the measurement matrix, by introduces additional
    % 1% of noise
    YM = kron(ym,ones(1,N)) + 0.01*kron(sigv,ones(1,N)).*randn(length(sigv),N);
    
    % True model propagation
    xtemp = A*xtrue(ii,:)' + B*u;
    xtrue(ii+1,:) = xtemp'; 
    
    % Estimated output for each ensemble
    for jj = 1:N
        estdopfreq = norm(X(3:4,jj))/c*f;
        r1e = X(1:2,jj) - towerloc(1,:)';
        r2e = X(1:2,jj) - towerloc(2,:)';
        Ye(:,jj) = [estdopfreq;doa(r1e);doa(r2e)];        
    end
    
    % Gain and Update
    dX = X*(eye(N)-1/N*ones(N));
    dY = Ye*(eye(N)-1/N*ones(N));    
    K = 1/(N-1)*dX*dY'*inv(1/(N-1)*dY*dY' + R);
    X = X+K*(YM-Ye);
    
    % Calculate the mean ensemble for output
    xe(ii,:) = mean(X,2);
    
    % Propagate all ensembles
    um = u+siga*randn(2,1);
    X = A*X + kron(B*um,ones(1,N));
    
    % --------------------------------------------------------------
    
    %% Extended Kalman filter
    % Generate measurement    
    dopfreqekf = norm(xekf(ii,3:4))/c*f;
    r1ekf = xekf(ii,1:2) - towerloc(1,:);
    r2ekf = xekf(ii,1:2) - towerloc(2,:);
    ye = [dopfreqekf;doa(r1ekf);doa(r2ekf)];  
    
    H = [0 0 xekf(ii,3)*f/(norm(xekf(ii,3:4))*c) xekf(ii,4)*f/(norm(xekf(ii,3:4))*c);
    -sin(doa(r1ekf))/norm(r1ekf) cos(doa(r1ekf))/norm(r1ekf) 0 0;
    -sin(doa(r2ekf))/norm(r2ekf) cos(doa(r2ekf))/norm(r2ekf) 0 0];
    
    K = Pekf*H'*inv(H*Pekf*H'+R);
    xekf(ii,:) = xekf(ii,:) + (K*(ym-ye))';
    Pekf = (eye(4)-K*H)*Pekf;
    
    % Propagation
    xekf(ii+1,:) = (A*xekf(ii,:)'+B*um)';
    Pefk = A*Pekf*A'+Qd;
end

figure(1)
plot(xe(:,1),xe(:,2))
hold on
plot(xtrue(:,1),xtrue(:,2),'r')
plot(xekf(:,1),xekf(:,2),'m')
xlabel('\fontsize{13}x-direction (km)')
ylabel('\fontsize{13}y-direction (km)')
grid on
set(gca,'fontsize',12)
plot(towerloc(1,1),towerloc(1,2),'x');
plot(towerloc(2,1),towerloc(2,2),'x');
