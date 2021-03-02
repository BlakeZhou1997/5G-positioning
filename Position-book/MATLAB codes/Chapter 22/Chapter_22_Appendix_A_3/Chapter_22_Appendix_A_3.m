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
gyro_y_filtered = filter(b_denoise,a_denoise,-gyro_X_Y_Z(:,2));
gyro_z_filtered = filter(b_denoise,a_denoise,-gyro_X_Y_Z(:,3));

acc_x_filtered = filter(b_denoise,a_denoise,acc_X_Y_Z(:,1));
acc_y_filtered = filter(b_denoise,a_denoise,-acc_X_Y_Z(:,2));
acc_z_filtered = filter(b_denoise,a_denoise,-acc_X_Y_Z(:,3));



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

Det_bias_acc = [AccxMean_bias;AccyMean_bias;AcczMean_bias];
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



%--------------------------------------------------------------------------
%      DETERMINISTIC SCALE-FACTOR REMOVAL
%-------------------------------------------------------------------------
%% Accelerometers deterministic Scale Factor
SF_AccX =   [0.9973 	0.9980      0.9988	 1.0000	 1.0010];
SF_AccY =   [0.9994 	1.0000      1.0008	 1.0019	 1.0030];
SF_AccZ =   [1.0002 	1.0008      1.0016   1.0020	 1.0025];

% Scale_factor variation due to the temperature
Temp = [-20 0 20 40 60];

sfAccX_pfit = polyfit(Temp, SF_AccX, 3);
sfAccY_pfit = polyfit(Temp, SF_AccY, 3);
sfAccZ_pfit = polyfit(Temp, SF_AccZ, 3);  

SF_Acc(1,:)     = polyval(sfAccX_pfit, temp(delay:delay+trans,1));
SF_Acc(2,:)     = polyval(sfAccY_pfit, temp(delay:delay+trans,1));
SF_Acc(3,:)     = polyval(sfAccZ_pfit, temp(delay:delay+trans,1));

% % 
SF_AccxMean = mean(SF_Acc(1,:));
SF_AccyMean = mean(SF_Acc(2,:));
SF_AcczMean =  mean(SF_Acc(3,:));
scale_factor_acc = ([SF_AccxMean;SF_AccyMean;SF_AcczMean]-ones(3,1));
% scale_factor_acc=[0;0;0];

%% Gryos deterministic Scale Factor
SF_GyrX = [	 1.0059	1.0042	1.0052	1.0033	1.0007];
SF_GyrY = [	 0.9988	0.9979	1.0001	0.9990	0.9980];
SF_GyrZ = [	 0.9993	1.0002	1.0028	1.0036	1.0042];

% Scale_factor variation due to the temperature
Temp = [-20 0 20 40 60];

sfGyroX_pfit = polyfit(Temp, SF_GyrX, 3);
sfGyroY_pfit = polyfit(Temp, SF_GyrY, 3);
sfGyroZ_pfit = polyfit(Temp, SF_GyrZ, 3);  
% 
SF_Gyro(1,:)     = polyval(sfGyroX_pfit, temp(delay:delay+trans,1));
SF_Gyro(2,:)     = polyval(sfGyroY_pfit, temp(delay:delay+trans,1));
SF_Gyro(3,:)     = polyval(sfGyroZ_pfit, temp(delay:delay+trans,1));

SF_GyroxMean = mean(SF_Gyro(1,:));
SF_GyroyMean = mean(SF_Gyro(2,:));
SF_GyrozMean =  mean(SF_Gyro(3,:));

scale_factor_gyros = [SF_GyroxMean;SF_GyroyMean;SF_GyrozMean]-ones(3,1);



%% INS NAVIGATION

%% GPS inizialization

load  '..\parser_log\GNSS_data.mat'; 
% GPS_pvt
% GPS_raw_data

    PosStart.Lat    = pvt(60).latitude;
    PosStart.Long   = pvt(60).longitude;
    PosStart.Alt    = pvt(60).height;
    PosStart.Vx     = pvt(60).v_lat;
    PosStart.Vy     = pvt(60).v_long;
    PosStart.Vz     = pvt(60).v_height;
        
    [ECEF]=llh2xyz([PosStart.Lat*pi/180;PosStart.Long*pi/180;PosStart.Alt]);
  
    PosStart.Xu =  ECEF(1);
    PosStart.Yu =  ECEF(2);
    PosStart.Zu =  ECEF(3);
    
%% INS Initialization
DURATION_TIME = 60; %s
n_samples = fs*DURATION_TIME;
PVT_rate = 10; %Hz
dt_medium_rate = 1/PVT_rate;

navdata.r_n          = zeros(3,floor(n_samples*(1/PVT_rate))); % Init O/P position LLH
navdata.r_e          = zeros(3,floor(n_samples*(1/PVT_rate))); % Init O/P position ECEF
navdata.v_n          = zeros(3,floor(n_samples*(1/PVT_rate))); % Init O/P velocity
navdata.v_e          = zeros(3,floor(n_samples*(1/PVT_rate))); % Init O/P velocity


%navdata.v_n(:,1)     = [v_N, v_E, v_D]'; % Init O/P velocity
navdata.r_n(:,1)     = [PosStart.Lat*pi/180,PosStart.Long*pi/180,PosStart.Alt]'; % Init O/P position
navdata.r_e(:,1)     = [PosStart.Xu;PosStart.Yu;PosStart.Zu];

navdata.R_n_e        = zeros(3,3,floor(n_samples*(1/PVT_rate)));
navdata.R_n_e(:,:,1) = pos2Cne(navdata.r_n(1,1),navdata.r_n(2,1));
navdata.v_e(:,1)     = [PosStart.Vx;PosStart.Vy;PosStart.Vz];
navdata.v_n(:,1)     = navdata.R_n_e(:,:,1).' *navdata.v_e(:,1);

navdata.Euler        = zeros(3,floor(n_samples*(1/PVT_rate))); % Init O/P attitude
navdata.Euler(3,1)   = heading;
navdata.Euler(2,1)   = pitch; % Init O/P attitude
navdata.Euler(1,1)   = roll; % Init O/P attitude

navdata.coning      = zeros(3,n_samples);
navdata.sculling    = zeros(3,n_samples);

navdata.R_b_n        = zeros(3,3,floor(n_samples*(1/PVT_rate)));
%navdata.R_b_n(:,:,1) = CBL;
navdata.R_b_n(:,:,1) = Euler2DCM(navdata.Euler(:,1));

navdata.R_b_e        = zeros(3,3,floor(n_samples*(1/PVT_rate))); % Init O/P DCM
navdata.R_b_e(:,:,1) = navdata.R_n_e(:,:,1)*navdata.R_b_n(:,:,1);

dt_INS = 1/fs;
iter_high_rate = 1;
iter_medium_rate = PVT_rate;
iteration_INS = 2;

angle_b_frame       = zeros(3,1);
delta_angle_b_frame = zeros(3,1);
vel_b_frame         = zeros(3,1);
delta_vel_b_frame   = zeros(3,1); 


for sample_idx=2:n_samples
         %% medium_rate part coning and sculling computation, Euler Angle computation, position and velocity estimation
             % coning and sculling
            [angle_b_frame,delta_angle_b_frame, navdata.coning(:,sample_idx), vel_b_frame,delta_vel_b_frame, navdata.sculling(:,sample_idx)]= coning_sculling_high_rate(...
                                                                                            [acc_x_filtered(delay+trans+sample_idx-1);acc_y_filtered(delay+trans+sample_idx-1);acc_z_filtered(delay+trans+sample_idx-1)].*g-Det_bias_acc,... 
                                                                                            [gyro_x_filtered(delay+trans+sample_idx-1);gyro_y_filtered(delay+trans+sample_idx-1);gyro_z_filtered(delay+trans+sample_idx-1)].*pi/180 - Det_bias_gyros*pi/180,vel_b_frame ,delta_vel_b_frame,...
                                                                                            angle_b_frame,delta_angle_b_frame,navdata.coning(:,sample_idx-1),... 
                                                                                            navdata.sculling(:,sample_idx-1), dt_INS,scale_factor_acc,scale_factor_gyros);

            iter_high_rate =   iter_high_rate + 1;                                                                  
             % computation of mechanization equation (attitude, pos and
             % velocity computation by the INS mechanization equations)
             
            if iter_high_rate == iter_medium_rate 
                        [navdata.R_b_e(:,:,iteration_INS), navdata.v_e(:,iteration_INS),navdata.r_e(:,iteration_INS) ]= ...
                                                                                strapdown_ecef_dcm_mediumRate2(navdata.R_b_e(:,:,iteration_INS-1),navdata.v_e(:,iteration_INS-1), navdata.r_e(:,iteration_INS-1),...
                                                                                                                vel_b_frame,navdata.sculling(:,sample_idx), dt_medium_rate,angle_b_frame,navdata.coning(:,sample_idx),navdata.r_n(:,iteration_INS-1));
                
                       navdata.r_n(:,iteration_INS)     = xyz2llh(navdata.r_e(:,iteration_INS));
                       navdata.R_n_e(:,:,iteration_INS) = pos2Cne(navdata.r_n(1,iteration_INS),navdata.r_n(2,iteration_INS));
                       navdata.v_n(:,iteration_INS)     = navdata.R_n_e(:,:,iteration_INS).'*navdata.v_e(:,iteration_INS);
           
                       navdata.R_b_n(:,:,iteration_INS) = navdata.R_n_e(:,:,iteration_INS).'*navdata.R_b_e(:,:,iteration_INS);                                                                                     
                       Euler = DCM2Euler(navdata.R_b_n(:,:,iteration_INS));      
                       navdata.Euler(:,iteration_INS) = [Euler.roll, Euler.pitch, Euler.yaw];                                                                                     
                                                                                                            
                       iteration_INS = iteration_INS+1;
                        
               
                        % reset for the high-rate part
                        
                         navdata.sculling(:, sample_idx) = zeros(3,1);
                         vel_b_frame = zeros(3,1);
                         delta_vel_b_frame = zeros(3,1);
                         iter_high_rate = 1;
                         angle_b_frame = zeros(3,1);
                         delta_angle_b_frame  = zeros(3,1);
                         navdata.coning(:,sample_idx)    = zeros(3,1);
            end
end

%% plot
figure,
subplot(311)
plot([1:length(navdata.r_e(1,1:end))],navdata.r_e(1,1:end))
hold on
plot([1:length(navdata.r_e(1,1:end))],repmat(PosStart.Xu,1,length(navdata.r_e(1,1:end))),'r')
title ('INS Solution Position ECEF')
xlabel('Num_samples')
ylabel ('ECEF-X [m]')
legend ('IMU', 'Reference')
subplot(312)
plot([1:length(navdata.r_e(2,1:end))],navdata.r_e(2,1:end))
hold on
plot([1:length(navdata.r_e(2,1:end))],repmat(PosStart.Yu,1,length(navdata.r_e(2,1:end))),'r')
xlabel('Num_samples')
ylabel ('ECEF-Y [m]')
legend ('IMU', 'Reference')
subplot(313)
plot([1:length(navdata.r_e(3,1:end))],navdata.r_e(3,1:end))
hold on
plot([1:length(navdata.r_e(3,1:end))],repmat(PosStart.Zu,1,length(navdata.r_e(3,1:end))),'r')
xlabel('Num_samples')
ylabel ('ECEF-Z [m]')
legend ('IMU', 'Reference')


figure,
subplot(211)
plot([1:length(navdata.r_e(1,1:end))],sqrt (navdata.r_e(1,1:end)-repmat(PosStart.Xu,1,length(navdata.r_e(1,1:end)))).^2+...
                                            (navdata.r_e(2,1:end)-repmat(PosStart.Yu,1,length(navdata.r_e(2,1:end)))).^2);
title ('INS Solution Horizontal Positional Error ECEF')
xlabel('Num_samples')
ylabel ('Horizontal Error [m]')
subplot(212)
plot([1:length(navdata.r_e(1,1:end))],sqrt (navdata.r_e(3,1:end)-repmat(PosStart.Zu,1,length(navdata.r_e(1,1:end)))).^2);
xlabel('Num_samples')
ylabel ('Vertical Error [m]')


figure,
subplot(311)
plot([1:length(navdata.v_e(1,1:end))],navdata.v_e(1,1:end))
hold on
plot([1:length(navdata.v_e(1,1:end))],repmat(PosStart.Vx,1,length(navdata.v_e(1,1:end))),'r')
title ('INS Solution Velocity ECEF')
xlabel('Num_samples')
ylabel ('ECEF-VelX [m/s]')
legend ('IMU', 'Reference')
subplot(312)
plot([1:length(navdata.v_e(2,1:end))],navdata.v_e(2,1:end))
hold on
plot([1:length(navdata.v_e(2,1:end))],repmat(PosStart.Vy,1,length(navdata.v_e(2,1:end))),'r')
xlabel('Num_samples')
ylabel ('ECEF-VelY [m/s]')
legend ('IMU', 'Reference')
subplot(313)
plot([1:length(navdata.v_e(3,1:end))],navdata.v_e(3,1:end))
hold on
plot([1:length(navdata.v_e(3,1:end))],repmat(PosStart.Vz,1,length(navdata.v_e(3,1:end))),'r')
xlabel('Num_samples')
ylabel ('ECEF-VelZ [m/s]')
legend ('IMU', 'Reference')



figure,
subplot(311)
plot([1:length(navdata.Euler(1,1:end))],navdata.Euler(1,1:end)*180/pi)
hold on
title ('INS Solution Euler')
xlabel('Num_samples')
ylabel ('Roll Angle [deg]')
subplot(312)
plot([1:length(navdata.Euler(1,1:end))],navdata.Euler(2,1:end)*180/pi)
xlabel('Num_samples')
ylabel ('Pitch Angle [deg]')
subplot(313)
plot([1:length(navdata.Euler(1,1:end))],navdata.Euler(3,1:end)*180/pi)
xlabel('Num_samples')
ylabel ('Yaw Angle [m/s]')

