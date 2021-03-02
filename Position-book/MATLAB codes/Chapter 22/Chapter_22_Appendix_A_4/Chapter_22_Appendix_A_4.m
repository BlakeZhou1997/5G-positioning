clear all
close all
clc

%'INS raw data are:
%acc_X_Y_Z [g]
%gyro_X_Y_Z [deg/sec]
%timestamp [sec]
%temp [deg]'

g = 9.80655; %[m/s^2];

load '..\parser_log\INS_data.mat'; 


%% ANTI_VIBRATION FILTER

fs = 100; % Hz
N=3; % order of the filter
Wn =5/fs; % cut-off frequency we want to cut everything up to 5Hz fi 
[b_denoise,a_denoise] = butter(N,Wn,'low'); 

gyro_x_filtered = filter(b_denoise,a_denoise,gyro_X_Y_Z(:,1));
gyro_y_filtered = filter(b_denoise,a_denoise,gyro_X_Y_Z(:,2));
gyro_z_filtered = filter(b_denoise,a_denoise,gyro_X_Y_Z(:,3));

acc_x_filtered = filter(b_denoise,a_denoise,acc_X_Y_Z(:,1));
acc_y_filtered = filter(b_denoise,a_denoise,acc_X_Y_Z(:,2));
acc_z_filtered = filter(b_denoise,a_denoise,acc_X_Y_Z(:,3));



%% Gyros deterministic Bias 
delay= 32;
TIME_CALIBRATION = 60; %sec
trans = fs*TIME_CALIBRATION;

gxMean = mean(gyro_x_filtered(delay:trans+delay,1));
gyMean = mean(gyro_y_filtered(delay:trans+delay,1));
gzMean = mean(gyro_z_filtered(delay:trans+delay,1));

Det_bias_gyros = [gxMean;gyMean;gzMean];

figure(1)
subplot(311)
plot(timestamp(32:end),gyro_x_filtered(32:end,1),'r');
hold on;
plot(timestamp(32:end),gyro_x_filtered(32:end,1)-Det_bias_gyros(1),'g');
ylabel('Gyro-X [deg/s]')
title('Gyroscopes data: RAW')
legend('filtered','bias-corrected')
subplot(312)
plot(timestamp(32:end),gyro_y_filtered(32:end,1),'r');
hold on;
plot(timestamp(32:end),gyro_y_filtered(32:end,1)-Det_bias_gyros(2),'g');
legend('filtered','bias-corrected')
ylabel('Gyro-Y [deg/s]')
subplot(313)
plot(timestamp(32:end),gyro_z_filtered(32:end,1),'r');
hold on;
plot(timestamp(32:end),gyro_z_filtered(32:end,1)-Det_bias_gyros(3),'g');
legend('filtered','bias-corrected')
ylabel('Gyro-Z [deg/s]')
xlabel('Time [s]')


%% Accelerometers deterministic Bias 

Bias_AccX = [0.0137     0.0038      -0.0070     -0.0181     -0.0268 ]*g;
Bias_AccY = [0.0137     0.0102      0.0044      4.46e-4 	-0.0033]*g;
Bias_AccZ = [0.0290 	5.9198e-04	-0.0248     -0.0312     -0.0674]*g;

Temp = [-20 0 20 40 60];

biasAccX_pfit = polyfit(Temp, Bias_AccX, 3);
biasAccY_pfit = polyfit(Temp, Bias_AccY, 3);
biasAccZ_pfit = polyfit(Temp, Bias_AccZ, 3);

AccX_bias   = polyval(biasAccX_pfit, temp(delay:delay+trans,1));
AccY_bias   = polyval(biasAccY_pfit, temp(delay:delay+trans,1));
AccZ_bias   = polyval(biasAccZ_pfit, temp(delay:delay+trans,1));

AccxMean_bias = mean(AccX_bias);
AccyMean_bias = mean(AccY_bias);
AcczMean_bias =  mean(AccZ_bias);

Det_bias_acc = [-AccxMean_bias;AccyMean_bias;-AcczMean_bias];
%% Attitude determination

g = 9.80655;  % gravitational constant  m/s^2
axMean = mean(acc_x_filtered(delay:trans+delay,1));
ayMean = mean(acc_x_filtered(delay:trans+delay,1));
azMean = mean(acc_x_filtered(delay:trans+delay,1));
AccMean = [axMean;ayMean;azMean]-Det_bias_acc;





pitch_two = asin(AccMean(1) / g);
roll_two = -asin(AccMean(2)/g*cos(pitch_two));
Cbl = zeros(3);
Cbl(1,1) = cos(pitch_two);
Cbl(1,2) = sin(pitch_two)*sin(roll_two);
Cbl(1,3) = sin(pitch_two)*cos(roll_two);
Cbl(2,2) = cos(roll_two);
Cbl(2,3) = -sin(roll_two);
Cbl(3,1) = -sin(pitch_two);
Cbl(3,2) =  cos(pitch_two)*sin(roll_two);
Cbl(3,3) = cos(pitch_two)*cos(roll_two);

CBL = Cbl;

[Angles]= DCM2Euler(CBL);
roll    = Angles.roll;
pitch   = Angles.pitch;
heading = Angles.yaw;