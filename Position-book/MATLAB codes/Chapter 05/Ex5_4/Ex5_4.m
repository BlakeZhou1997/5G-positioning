clear all
close all
clc

% Simulate and estimate the endoscopy capsule position in human body
% Compare both UKF and EKF
% Required matlab file
%   endoscopy_position.m, gen_dbs_pos.m


% Constat Parameters
d2r = pi/180;
c = 300*10^8; % light speed cm/sec
m = 68;
dt = 1;
tp = 0.001;
Re = [200 100 50]; % position of Emitter which is fixukfd.


%% Generate 16 receiver (or dbs) locations
Rd = gen_dbs_pos(4,4,40);
dbs_no = size(Rd,1); % number of DBS

% Simulate Endoscopy position, velocity and acceleration
endoscopy_position;
v(1,:)=[0 0 0];
for k=1:m-1
    a(k,:)=2*(R(k+1,:)-R(k,:) - v(k,:)*dt)/dt^2;
    v(k+1,:)=a(k,:)*dt+v(k,:);
end


% Dynamis system model
A=[1 0 0 dt 0 0;
    0 1 0 0 dt 0;
    0 0 1 0 0 dt;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
B=[0.5*dt^2 0 0;
    0   0.5*dt^2 0;
    0   0    0.5*dt^2;
    dt   0    0;
    0     dt   0;
    0      0    dt];


% Measurement Noise Setup
sigt = 0.1; %TOA error
sigr = 0.1*d2r; %DOA error
errmat = kron([sigt sigr sigr]',ones(1,dbs_no));
% Covariance Matrices
RR = diag(reshape(errmat,1,3*dbs_no)); %measurement noise covariance alpha, beta and gamma are the proposed measurement noises
% Process Noise Covariance
q = 0.01*9.81;
Q = q^2*[(1/3)*(dt^3)*eye(3) (1/2)*(dt^2)*eye(3); (1/2)*(dt^2)*eye(3) dt*eye(3)];




for imonte = 1:20
    % initialization for EKF state vector and Covariance
    xekf = [0 8 0 0 0 0];
    Pekf = [5 5 5  .1 .1 .1]; %state covariance
    pekf = diag(Pekf);
    
    % initialization for UKF state vector and Covariance
    xukf = [0 8 0 0 0 0];
    Pukf = [5 5 5  .1 .1 .1]; %state covariance
    pukf = diag(Pukf);
    
    % unscented parameters
    n = length(xukf);
    kappa = 3 - n;
    alpha = 0.1;
    beta = 0.5;
    lam = alpha^2*(n+kappa) - n;
    
    w0m = lam/(n+ lam);
    wim = 1/(2*(n+ lam));
    
    w0c = lam/(n+lam)+(1-alpha^2+beta);
    wic = 1/(2*(n+lam));
    
    for k=1:size(R,1)-1
        %measurement with noise
        errm = errmat'.*randn(dbs_no,3); %output noise data
        anoise(k,:) = q*randn(1,3);
        u=a(k,:)'+anoise(k,:)'; % accelerometer noise data
        ym = [];
        gam0 = [];
        
        
        sig = chol((n + lam)* (pukf + Q));
        xx0 = xukf(k,:)';
        xx = kron(xx0,ones(1,2*n)) + [sig -sig];
        
        % Propagation for xx0
        xx0 = A*xx0+B*u;
        for jj=1: dbs_no
            yy =[norm(Re-R(k+1,:))-norm(R(k+1,:)-Rd(jj,:))+tp*c;atan2(R(k+1,2)-Rd(jj,2),R(k+1,1)-Rd(jj,1));atan2(R(k+1,3)-Rd(jj,3),sqrt((R(k+1,1)-Rd(jj,1))^2+(R(k+1,2)-Rd(jj,2))^2))]+errm(jj,:)';
            ym = [ym;yy];
            %% with gyro
            % estimated measurement xx0
            yy =[(norm(Re-xx0(1:3)')-norm(xx0(1:3)'-Rd(jj,:)))+tp*c;atan2(xx0(2)-Rd(jj,2),xx0(1)-Rd(jj,1));atan2(xx0(3)-Rd(jj,3),sqrt((xx0(1)-Rd(jj,1))^2+(xx0(2)-Rd(jj,2))^2))];
            gam0 = [gam0;yy];
        end
        
        
        for ii = 1:2*n
            % Propagation for xx
            xx(:,ii) = A*xx(:,ii)+B*u;
            % Estimated measurement for xx
            gamtemp = [];
            for jj=1: dbs_no
                yy =[(norm(Re-xx(1:3,ii)')-norm(xx(1:3,ii)'-Rd(jj,:)))+tp*c;atan2(xx(2,ii)-Rd(jj,2),xx(1,ii)-Rd(jj,1));atan2(xx(3,ii)-Rd(jj,3),sqrt((xx(1,ii)-Rd(jj,1))^2+(xx(2,ii)-Rd(jj,2))^2))];
                gamtemp = [gamtemp;yy];
            end
            gamma(:,ii) = gamtemp;
        end
        
        % calculate mean
        xukf(k+1,:)=w0m*xx0'+wim*sum(xx,2)';
        yukf=w0m*gam0+wim*sum(gamma,2);
        
        % Calculate covariances
        ppx0=w0c*(xx0'-xukf(k+1,:))'*(xx0'-xukf(k+1,:));
        ppmatx=xx-kron(xukf(k+1,:)',ones(1,2*n));
        pukf=ppx0+wic*ppmatx*ppmatx'+Q;
        
        ppy0=w0c*(gam0-yukf)*(gam0-yukf)';
        ppmaty=gamma-kron(yukf,ones(1,2*n));
        pyy=ppy0+wic*ppmaty*ppmaty';
        pvv=pyy+RR;
        
        pxy=w0c*(xx0'-xukf(k+1,:))'*(gam0-yukf)'+wic*ppmatx*ppmaty';
        
        
        % Kalman Gain Calculation
        K = pxy*inv(pvv);
        %------------------------------------------
        % Update
        %------------------------------------------
        pukf=pukf-K*pvv*K';
        Pukf(k+1,:)=diag(pukf);
        
        dy = ym-yukf;
        
        dy = reshape(dy,3,dbs_no);
        dyr = dy(2:3,:);
        ind = find(abs(dyr) > pi);
        if length(ind) > 0
            dyr(ind) = dyr(ind)-sign(dyr(ind))*2*pi;
            dy(2:3,:) = dyr;
        end
        dy(2:3,:) = dyr;
        dy = reshape(dy,dbs_no*3,1);
        xukf(k+1,:)=xukf(k+1,:)+(K*dy)';
        
        
        %% ------------------------------------------
        % Extended Kalman Filter
        %------------------------------------------
        %------------------------------------------
        %Estimated State Propagation
        %------------------------------------------        
        xekf(k+1,:) = (A*xekf(k,:)'+B*u)';
        pekf = A*pekf*A'+Q;
        Pekf(k+1,:)=diag(pekf);
        
        yekf = [];
        Hekf = [];
        for jj=1: dbs_no
            
            % estimated measurement
            yy =[(norm(Re-xekf(k+1,1:3))-norm(xekf(k+1,1:3)-Rd(jj,:)))+tp*c;atan2(xekf(k+1,2)-Rd(jj,2),xekf(k+1,1)-Rd(jj,1));atan2(xekf(k+1,3)-Rd(jj,3),sqrt((xekf(k+1,1)-Rd(jj,1))^2+(xekf(k+1,2)-Rd(jj,2))^2))];
            yekf = [yekf;yy];
            %------------------------------------------
            % H, Sensitivity Matrix Calculation
            %------------------------------------------
            b1=norm(Re-xekf(k+1,1:3));
            b2=norm(xekf(k+1,1:3)-Rd(jj,:));
            a1=(-1/b1)*(Re(1)-xekf(k+1,1))-(-1/b2)*(Rd(jj,1)-xekf(k+1,1));
            a2=(-1/b1)*(Re(2)-xekf(k+1,2))-(-1/b2)*(Rd(jj,2)-xekf(k+1,2));
            a3=(-1/b1)*(Re(3)-xekf(k+1,3))-(-1/b2)*(Rd(jj,3)-xekf(k+1,3));
            c1=-sin(yy(2))/cos(yy(3))/b2;
            c2=cos(yy(2))/cos(yy(3))/b2;
            c3=0;
            d1=-cos(yy(2))*sin(yy(3))/b2;
            d2=-sin(yy(2))*sin(yy(3))/b2;
            d3=cos(yy(3))/b2;
            hh = [a1,a2,a3,0,0,0;
                c1,c2,c3,0,0,0;
                d1,d2,d3,0,0,0];
            Hekf = [Hekf;hh];
        end
        
        %Kalman Filter gain,updates and propagation
        %------------------------------------------
        % Kalman Gain Calculation
        %------------------------------------------
        gain = pekf*Hekf'*inv(Hekf*pekf*Hekf'+RR);
        %------------------------------------------
        % Update
        %------------------------------------------
        pekf = (eye(6)-gain*Hekf)*pekf;
        dy = ym-yekf;
        dy = reshape(dy,3,dbs_no);
        dyr = dy(2:3,:);
        ind = find(abs(dyr) > pi);
        if length(ind) > 0
            dyr(ind) = dyr(ind)-sign(dyr(ind))*2*pi;
            dy(2:3,:) = dyr;
        end
        dy(2:3,:) = dyr;
        dy = reshape(dy,dbs_no*3,1);
        xekf(k+1,:)=xekf(k+1,:)+((gain*(dy)))';
       
        
    end
    
    errukf = xukf(2:end,1:3)-R(2:end,:);
    errekf = xekf(2:end,1:3)-R(2:end,:);
    ukfrmse(imonte,:) = sum(errukf.^2,2).^.5/sqrt(3);
    ekfrmse(imonte,:) = sum(errekf.^2,2).^.5/sqrt(3);
end

meanukfrmse = mean(ukfrmse,1);
meanekfrmse = mean(ekfrmse,1);
% Three sigma
sig3 = 3*Pukf.^.5;
% create a x-vector for plotting
t=[1:size(R,1)];
figure(1)
subplot(311)
plot(t,(xukf(:,1)-R(:,1)),t,sig3(:,1)/4,t,-sig3(:,1)/4)
ylabel('\fontsize{14}x-err.(cm)')
xlabel('\fontsize{14} Time(min)')
set(gca,'fontsize',12)
grid on
xlim([1 68])
ylim([-0.5 0.5])
subplot(312)
plot(t,(xukf(:,2)-R(:,2)),t,sig3(:,2)/4,t,-sig3(:,2)/4)
ylabel('\fontsize{14}y-err.(cm)')
xlabel('\fontsize{14}Time (min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])
subplot(313)
plot(t,(xukf(:,3)-R(:,3)),t,sig3(:,3)/4,t,-sig3(:,3)/4)
ylabel('\fontsize{12}z-rr.(cm)')
xlabel('\fontsize{12}Time(min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])

% plot EKF 3-axis error (option)
figure(2)
ekfsig3 = 3*Pekf.^.5;
subplot(311)
plot(t,(xekf(:,1)-R(:,1)),t,ekfsig3(:,1)/4,t,-ekfsig3(:,1)/4)
ylabel('\fontsize{14}x-err.(cm)')
xlabel('\fontsize{14} Time(min)')
set(gca,'fontsize',12)
grid on
xlim([1 68])
ylim([-0.5 0.5])
subplot(312)
plot(t,(xekf(:,2)-R(:,2)),t,ekfsig3(:,2)/4,t,-ekfsig3(:,2)/4)
ylabel('\fontsize{14}y-err.(cm)')
xlabel('\fontsize{14}Time (min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])
subplot(313)
plot(t,(xekf(:,3)-R(:,3)),t,ekfsig3(:,3)/4,t,-ekfsig3(:,3)/4)
ylabel('\fontsize{12}z-rr.(cm)')
xlabel('\fontsize{12}Time(min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])

%Plot RMSE comparison
figure(3)
plot(t(2:end),meanekfrmse,t(2:end),meanukfrmse)
legend('EKF','UKF')
xlabel('\fontsize{13}Samples')
ylabel('\fontsize{13}Mean RMSE (cm)')
grid on
set(gca,'fontsize',12)
xlim([1 68])