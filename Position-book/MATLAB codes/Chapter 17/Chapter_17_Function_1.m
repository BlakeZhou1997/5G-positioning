function [Rx_pos,true_toa,meas_toa]=...
    Chapter_17_Function_1(Ktheta,noise_var,epsilon,kappa)
% INPUT:
%   Ktheta      : no of generated TOA measurements
%   noise_var   : variance of measurement noise [in microsec]
%   epsilon     : amount of contamination
%   kappa       : parameter for epsilon-contaminated Gaussian mixture 
% OUTPUT:
%   true_toa    : actual TOA
%   meas_toa    : noisy TOA

c = 3e8;       % speed of light [m/s^2]
Rx_x0 = 10*1e3; Rx_y0 = 10*1e3; % initial position of sensor [m]
Tx_x0 = 5*1e3; Tx_y0 = 5*1e3; % initial position of MT [m]
theta_step = 5; 
% calculate moving positions
Rx_x1 = []; Rx_y1 = []; count = 1;
for theta = 45:theta_step:290;
    Rx_x1(count) = sqrt(Rx_x0^2 + Rx_y0^2)*cos(theta*pi/180);
    Rx_y1(count) = sqrt(Rx_x0^2 + Rx_y0^2)*sin(theta*pi/180);
    count=count+1;
end
Rx_pos = [Rx_x1.' Rx_y1.'];
Range = sqrt((Tx_x0-Rx_x1).^2 + (Tx_y0-Rx_y1).^2); 
true_toa = Range.'/c; % in sec

thermal_noise = randn(Ktheta,1)*noise_var*1e-6; % meas error [sec]

%% generate NLOS errors from epsilon-contaminated Gaussian mixture
%two gaussian distributions are built, one with sigma1; one with sigma2
x1 = randn(Ktheta,1);
x2 = sqrt(kappa).*randn(Ktheta,1);
% creates random number between 0 and 1 of a uniform distribution 
mixmat = rand(Ktheta,1);
%creates a mask, mixsel1,2 takes values of zero or one
mixsel1 = mixmat>epsilon;
mixsel2 = mixmat<epsilon;
%random gaussian numbers are multiplied with 1 or 0 in order to obtain 
%the overall distribution 
eCMM = abs(mixsel1.*x1 + mixsel2.*x2); % in microsecond
eCMM = eCMM.*1e-6; % in sec

%% random position for the NLOS errors
binmat = zeros(Ktheta,1); impul_sam = randperm(Ktheta-4)+2;
for count = 1:5
    idx = impul_sam(count);
    binmat(idx-2:idx+2) = 1;
end

meas_toa = true_toa + thermal_noise + eCMM.*binmat; % dist_err.*binmat; %