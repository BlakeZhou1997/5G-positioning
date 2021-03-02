clear all
close all
clc

% Simulate and estimate the endoscopy capsule position in human body
%
% Required matlab file
%   endoscopy_position.m, gen_dbs_pos.m


% Constat Parameters
d2r = pi/180;
c = 300*10^8; % light speed cm/sec
m = 68;
dt = 1;
tp = 0.001;
Re = [200 100 50]; % position of Emitter which is fixed.


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
anoise = q*randn(m,3);
Q = q^2*[(1/3)*(dt^3)*eye(3) (1/2)*(dt^2)*eye(3); (1/2)*(dt^2)*eye(3) dt*eye(3)];


% initialization for EKF state vector and Covariance
xe = [0 8 0 0 0 0];
P = [5 5 5  .1 .1 .1]; %state covariance
pcov = diag(P);

for k=1:size(R,1)-1
    %measurement with noise
    errm = errmat'.*randn(dbs_no,3);
    %         errm = 0.*randn(dbs_no,3);
    ym = [];
    ye = [];
    H = [];
    
    for jj=1: dbs_no
        
        yy =[norm(Re-R(k,:))-norm(R(k,:)-Rd(jj,:))+tp*c;atan2(R(k,2)-Rd(jj,2),R(k,1)-Rd(jj,1));atan2(R(k,3)-Rd(jj,3),sqrt((R(k,1)-Rd(jj,1))^2+(R(k,2)-Rd(jj,2))^2))]+errm(jj,:)';
        ym = [ym;yy];
        
        %% with gyro
        % estimated measurement
        yy =[(norm(Re-xe(k,1:3))-norm(xe(k,1:3)-Rd(jj,:)))+tp*c;atan2(xe(k,2)-Rd(jj,2),xe(k,1)-Rd(jj,1));atan2(xe(k,3)-Rd(jj,3),sqrt((xe(k,1)-Rd(jj,1))^2+(xe(k,2)-Rd(jj,2))^2))];
        ye = [ye;yy];
        %------------------------------------------
        % H, Sensitivity Matrix Calculation
        %------------------------------------------
        b1=norm(Re-xe(k,1:3));
        b2=norm(xe(k,1:3)-Rd(jj,:));
        a1=(-1/b1)*(Re(1)-xe(k,1))-(-1/b2)*(Rd(jj,1)-xe(k,1));
        a2=(-1/b1)*(Re(2)-xe(k,2))-(-1/b2)*(Rd(jj,2)-xe(k,2));
        a3=(-1/b1)*(Re(3)-xe(k,3))-(-1/b2)*(Rd(jj,3)-xe(k,3));
        c1=-sin(yy(2))/cos(yy(3))/b2;
        c2=cos(yy(2))/cos(yy(3))/b2;         
        c3=0;
        d1=-cos(yy(2))*sin(yy(3))/b2;
        d2=-sin(yy(2))*sin(yy(3))/b2;
        d3=cos(yy(3))/b2;
        hh = [a1,a2,a3,0,0,0;
            c1,c2,c3,0,0,0;
            d1,d2,d3,0,0,0];
        H = [H;hh];
    end
    
    %Kalman Filter gain,updates and propagation
    %------------------------------------------
    % Kalman Gain Calculation
    %------------------------------------------
    gain = pcov*H'*inv(H*pcov*H'+RR);
    %------------------------------------------
    % Update
    %------------------------------------------
    pcov = (eye(6)-gain*H)*pcov;
    dy = ym-ye;
    dy = reshape(dy,3,dbs_no);
    dyr = dy(2:3,:);
    ind = find(abs(dyr) > pi);
    if length(ind) > 0
        dyr(ind) = dyr(ind)-sign(dyr(ind))*2*pi;
        dy(2:3,:) = dyr;
    end
    dy(2:3,:) = dyr;
    dy = reshape(dy,dbs_no*3,1);
    xe(k,:)=xe(k,:)+((gain*(dy)))';
    %------------------------------------------
    %Estimated State Propagation
    %------------------------------------------
    a(k,:)=a(k,:)+anoise(k,:);
    xe(k+1,:) = (A*xe(k,:)'+B*(a(k,:)'))';
    pcov = A*pcov*A'+Q;
    P(k+1,:)=diag(pcov);
    
end


% Three sigma
sig3 = 3*P.^.5;
%Error
err = xe(:,1:3) -R;
% create a x-vector for plotting
t=[1:size(R,1)];
figure(1)
subplot(311)
plot(t,(xe(:,1)-R(:,1)),t,sig3(:,1)/4,t,-sig3(:,1)/4)
ylabel('\fontsize{14}x-err.(cm)')
xlabel('\fontsize{14} Time(min)')
set(gca,'fontsize',12)
grid on
xlim([1 68])
ylim([-0.5 0.5])
subplot(312)
plot(t,(xe(:,2)-R(:,2)),t,sig3(:,2)/4,t,-sig3(:,2)/4)
ylabel('\fontsize{14}y-err.(cm)')
xlabel('\fontsize{14}Time (min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])
subplot(313)
plot(t,(xe(:,3)-R(:,3)),t,sig3(:,3)/4,t,-sig3(:,3)/4)
ylabel('\fontsize{12}z-rr.(cm)')
xlabel('\fontsize{12}Time(min)')
set(gca,'fontsize',12)
xlim([1 68])
ylim([-0.5 0.5])



