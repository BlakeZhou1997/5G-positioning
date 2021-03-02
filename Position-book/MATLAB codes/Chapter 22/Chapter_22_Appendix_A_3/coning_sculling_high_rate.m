function [alfa,delta_alfa, coning, vel,delta_vel, sculling]=coning_sculling_high_rate(a,w,vel_prec,delta_vel_prec,...
                                                           alfa_prec,delta_alfa_prec,coning_prec, sculling_prec, dt_high_rate,sa,sg)

 %% INPUT
 % a                    --> DE-NOISED Accelerometers raw data array (3x1) corrected with bias [m/s^2] at time n time high_rate (100Hz)
 % w                    --> DE-NOISED Angular rate array (3x1) corrected with gyro bias [rad/s] at time  n high_rate (100Hz) 
 %vel_prec              --> vel_b_frame at time n-1 high_rate [m/s]
 %vel_prec              --> delta_vel_b_frame at time n-1 high_rate [m/s] 
 % alfa_prec            --> angle_b_frame at time n-1  high_rate [rad]
 % delta_alfa_prec      --> delta angle_b_frame at time n-1 high_rate [rad]
 % coning_prec          --> coning at time n-1 high_rate
 % sculling_prec        --> sculling at time n-1 high_rate
 % dt_high_rate         --> dt_INS [100Hz]
 %sa                    --> scale_factor accelerometers array(3x1)
 %sg                    --> scale_factor gyroscopes array(3x1)
  
% only the sculling and coning is computed
%% coning increment every 100Hz
% increment due to the rotation
delta_alfa = w*dt_high_rate;
delta_alfa =([1/(1.00+sg(1)) 0 0; 0 1/(1.00+sg(2)) 0; 0 0 1/(1.00+sg(3))])*delta_alfa;
alfa =  alfa_prec+delta_alfa;
% increment of coning SAVAGE paper eq 45
delta_coning = cross(0.5*(alfa_prec+(1/6)*delta_alfa_prec),delta_alfa);
coning = coning_prec+delta_coning;


%% sculling increment every 100Hz
% increment velocity
delta_vel = a*dt_high_rate;
delta_vel =([1/(1.00+sa(1)) 0 0; 0 1/(1.00+sa(2)) 0; 0 0 1/(1.00+sa(3))])*delta_vel;
vel = vel_prec+delta_vel;
%increment of sculling SAVAGE paper eq 61
i = cross((alfa_prec+(1/6)*delta_alfa_prec),delta_vel);
ii = cross((vel_prec+(1/6)*delta_vel_prec),delta_alfa);
delta_sculling = 0.5*(i+ii);
sculling = sculling_prec+delta_sculling;

