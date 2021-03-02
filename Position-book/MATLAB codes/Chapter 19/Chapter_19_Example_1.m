%Matlab code for example 19.1:

%initialization for states and time
T = 2*3600; %orbit period, 2hrs
n = 2*pi/T; %mean motion
 
%initial position and velocity
r = [200 300]; 
v = [.015 -2*n*r(1)]; %velocity is required to fulfill the cartesian initial condition
x = [r v];
 
%Total simulation period
dt = 10;
tf = 6*3600;
t = 0:dt:tf;
m = length(t);
 
%Relative attitude, assume same for all the time
C = [cos(45*pi/180) sin(45*pi/180);
     -sin(45*pi/180) cos(45*pi/180)];
 
%Generate truth position
for i = 2:m
    [tt,xx]=ode45(@CW2D_ODE,[t(i-1) t(i)],x(i-1,:),[],n);
    x(i,:) = xx(end,:);
end
 
%Measurement errors
sigr = 1;
sigp = 17e-5;
errm = [sigr*randn(m,2)];
 
%Estimated states, states covariance, noise covariance initialization
xe = x(1,:); %estimated states
P = [5 5 1 1]; %state covariance
pcov = diag(P);
R = diag([sigr^2 sigr^2]); %measurement noise covariance
Q = diag([sqrt(10)*10^(-11)*ones(1,2)]).^2; %process noise covariance
G = [zeros(2);eye(2)];
F = [0 0 1 0;
    0 0 0 1;
    3*n^2 0 0 2*n;
     0 0 -2*n 0 ];
 
for i=1:m-1
    %measurement with noise
    ym = [x(i,1:2)']+errm(i,:)';
    
    %Kalman Filter algorithm
    %estimated measurement
    ye = [xe(i,1:2)'];
    
    
    %Kalman Filter gain and updates
    H = [eye(2) zeros(2)];
    gain = pcov*H'*inv(H*pcov*H'+R);
    pcov = (eye(4)-gain*H)*pcov;
    xe(i,:)=xe(i,:)+(gain*(ym-ye))';
    
    %------------------------------------------
    %Estimated State Propagation
    %------------------------------------------
    [tt,xx]=ode45(@CW2D_ODE,[t(i) t(i+1)],xe(i,:),[],n);
    xe(i+1,:) = xx(end,:);
    
    P(i,:) = diag(pcov);
    
    covmat = [F G [Q;zeros(2)]];
    pmat = reshape(pcov,1,16); %reshape 4x4 matrix into 1x16 vector
    [tt,xx] = ode45(@covprop,[t(i) t(i+1)],pmat,[],covmat);
    pcov = reshape(xx(end,:),4,4); %reshape 1x16 vector into 4x4 matrix
    P(i+1,:) = diag(pcov);
end
 
sig3 = 3*P.^.5;
 
absrms = (sum((xe(:,1:2)-x(:,1:2)).^2,2)).^.5;    
figure(1)
subplot(211)
plot(t/60,xe(:,1)-x(:,1),t/60,sig3(:,1),t/60,-sig3(:,1))
ylabel('\fontsize{14}X-position error (m)')
xlabel('\fontsize{14}Time (min)')
axis([0 360 -1 1])
subplot(212)
plot(t/60,xe(:,2)-x(:,2),t/60,sig3(:,2),t/60,-sig3(:,2))
ylabel('\fontsize{14}Y-position error (m)')
xlabel('\fontsize{14}Time (min)')
axis([0 360 -1 1])
 
figure(2)
subplot(211)
plot(t/60,xe(:,3)-x(:,3),t/60,sig3(:,3),t/60,-sig3(:,3))
ylabel('\fontsize{14}X-velocity error (m/s)')
xlabel('\fontsize{14}Time (min)')
axis([0 360 -.001 .001])
subplot(212)
plot(t/60,xe(:,4)-x(:,4),t/60,sig3(:,4),t/60,-sig3(:,4))    
ylabel('\fontsize{14}Y-velocity error (m/s)')
xlabel('\fontsize{14}Time (min)')
axis([0 360 -.001 .001])

%%%%%%%%%%%%%%%%
function f=CW2D_ODE(t,x,n)
 
f=[0 0 0 0]';
 
f(1)=x(3);
f(2)=x(4);
f(3)= 2*n*x(4)+3*n^2*x(1);
f(4)=-2*n*x(3);
