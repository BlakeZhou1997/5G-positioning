clear all
close all
clc


% Required matlab files and data files
%   almyuma_520_dat.alm
%   read_yuma.m, alm2satposvel.m, orb2eci.m

% Constant Parameters
global CONST_MU_E CONST_OMEGA_E
d2r = pi/180; % deg to radian
mu = 398600.4; % earth gravitational parameter
CONST_MU_E = mu; % required by alm2satposvel.m
CONST_OMEGA_E = 7292115.1467e-11; % earth rotation rate

% Time configuration
dt = 1; %sampling time
tf = 6000;
t = 0:dt:tf;


% Read GPS satellite position
gps_time = (8+31*365)*3600*24;
alm_param=read_yuma('almyuma_520_dat.alm');
gpsnum = 4;

% noise parameter
siggps = 0.1; % gps range error std in km
sigq = 0.001; % process noise std, in km
R = siggps^2*eye(gpsnum); % standard KF noise covariance, we always use 4 gps
Rckf = [siggps^2*eye(gpsnum) zeros(gpsnum,3); zeros(3,gpsnum) 1e-6*eye(3)];
Q = sigq^2*eye(3);
D = [zeros(3);eye(3)];


% Repeat 40 monte carlo simulation
for imonte = 1:20
    % True satellite
    oev = [6978.145 0.0 98*d2r 0*d2r 100*d2r 15*d2r];
    [r,v] = orb2eci(oev);
    xtrue = [r' v'];
    dconst = cpm(r)*v/norm(cpm(r)*v); % constraint condition
    
    % CKF initialization
    oev = [6978.145 0.0 97.8*d2r 0*d2r 100*d2r 15.2*d2r];
    [r,v] = orb2eci(oev);
    xckf = [r' v'];
    pcov = [100 100 100 1 1 1];
    Pckf = diag(pcov).^2;
    
    % EKF initialization
    xekf = xckf;
    Pekf = Pckf;
    pekf = pcov;
    for ii = 1:dt:length(t)-1
        % Generate GPS satellite position
        [prn,sv_xyz,sv_xyz_dot]=alm2satposvel(gps_time+t(ii), alm_param);
        
        allrange = sum((sv_xyz - kron(xtrue(ii,1:3),ones(length(prn),1))).^2,2).^.5;
        datamat = [allrange (1:length(prn))']; % do not use prn vector because some gps may not available
        tempdata = sortrows(datamat,1); % sort the matrix from closest range to further
        
        ym = tempdata(1:gpsnum,1) + siggps*randn(gpsnum,1);
        gpsprn = tempdata(1:gpsnum,2);
        gpspos = sv_xyz(gpsprn,:);
        
        %-------------------------------------------------------------------
        %CKF Process
        %-------------------------------------------------------------------
        % CKF Update
        ye = sum((gpspos - kron(xckf(ii,1:3),ones(gpsnum,1))).^2,2).^.5; % GPS measurement
        C = [(kron(xckf(ii,1:3),ones(gpsnum,1))-gpspos)./kron(ye,ones(1,3)) zeros(gpsnum,3)];
        
        % Constrain update
        estb = cpm(xckf(ii,1:3))*xckf(ii,4:6)'/norm(cpm(xckf(ii,1:3))*xckf(ii,4:6)'); %Constraint Estimation
        r = xckf(ii,1:3)'; v = xckf(ii,4:6)';
        dbdx = [-cpm(v)/norm(cpm(r)*v) + 0.5*cpm(r)*v*(cpm(cpm(r)*v)*v)'/norm(cpm(r)*v)^(3/2) cpm(r)/norm(cpm(r)*v)-0.5*cpm(r)*v*(cpm(cpm(r)*v)*r)'/norm(cpm(r)*v)^(3/2)];
        
        ymckf = [ym;dconst];
        yeckf = [ye;estb];
        Cckf = [C;dbdx];
        
        K = Pckf*Cckf'*inv(Cckf*Pckf*Cckf' + Rckf);
        xckf(ii,:) = xckf(ii,:) + (K*(ymckf - yeckf))';
        Pckf = (eye(6) - K*Cckf)*Pckf;
        pcov(ii,:) = diag(Pckf);
        
        
        % Propagation
        x0 = xckf(ii,:)';
        [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
        xckf(ii+1,:) = xx(end,:);
        
        F = [zeros(3) eye(3);
            diff_orbit(xckf(ii,1:3),mu) zeros(3);];
        p0 = reshape(Pckf,6*6,1);
        [~,xx] = ode45(@(tt,pl)covprop(tt,pl,F,D,Q),[0 dt],p0);
        pfinal = xx(end,:);
        Pckf = reshape(pfinal,6,6);
        pcov(ii+1,:) = diag(Pckf);
        
        %-------------------------------------------------------------------
        %EKF Process
        %-------------------------------------------------------------------
        % EKF Update
        ye = sum((gpspos - kron(xekf(ii,1:3),ones(gpsnum,1))).^2,2).^.5;
        C = [(kron(xekf(ii,1:3),ones(gpsnum,1))-gpspos)./kron(ye,ones(1,3)) zeros(gpsnum,3)];
        K = Pekf*C'*inv(C*Pekf*C' + R);
        xekf(ii,:) = xekf(ii,:) + (K*(ym - ye))';
        Pekf = (eye(6) - K*C)*Pekf;
        pekf(ii,:) = diag(Pekf);
        
        x0 = xekf(ii,:)';
        [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
        xekf(ii+1,:) = xx(end,:);
        
        F = [zeros(3) eye(3);
            diff_orbit(xekf(ii,1:3),mu) zeros(3);];
        p0 = reshape(Pekf,6*6,1);
        [~,xx] = ode45(@(tt,pl)covprop(tt,pl,F,D,Q),[0 dt],p0);
        pfinal = xx(end,:);
        Pekf = reshape(pfinal,6,6);
        pekf(ii+1,:) = diag(Pekf);
        
        %-------------------------------------------------------------------
        %Truth States Propagation
        %-------------------------------------------------------------------
        x0 = xtrue(ii,:)';
        [~,xx] = ode45(@(tt,x)orbit(tt,x,mu),[0 dt],x0);
        xtrue(ii+1,:) = xx(end,:);
    end
    ckferr = xckf - xtrue;
    ekferr = xekf - xtrue;
    ckfrmsemonte(:,imonte) = sum(ckferr(:,1:3).^2,2).^.5/sqrt(3);
    ekfrmsemonte(:,imonte) = sum(ekferr(:,1:3).^2,2).^.5/sqrt(3);
    
end

ckfrmse = mean(ckfrmsemonte,2);
ekfrmse = mean(ekfrmsemonte,2);

figure(1)
plot(t/60,ckfrmse*1000,t/60,ekfrmse*1000)
grid on
xlabel('\fontsize{13}Time (min)')
ylabel('\fontsize{13}RMSE (m)')
legend('CKF','EKF')
set(gca,'fontsize',12)