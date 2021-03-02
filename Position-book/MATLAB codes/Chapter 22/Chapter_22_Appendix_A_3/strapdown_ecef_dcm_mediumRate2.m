function [Cbe_new, Ve_new, ecef_new]=strapdown_ecef_dcm_mediumRate2(Cbe, Ve,ecef,...
                                                                        delta_vel,sculling, dt,rotation_angle,coning,pos_ned)

 %% INPUT
%  Cbe                 --->  Attitude DCM Matrix (3x3) from body to ECEF frame at time medium_rate -1
%  Ve                  --->  Velocity Array (3x1) in ECEF frame at time medium_rate -1
%  ecef                --->  Position Array (3x1) in ECEF frame at time mediumu_rate -1
%  delta_vel           --->  vel_b_frame at time n high_rate [m/s]
%  sculling            --->  sculling at time n high_rate
%  dt                  --->  integration time at medium rate (e.g. 10 Hz) [s]
% rotation_angle       --->  angle_b_frame at time n high_rate [rad]
%  coning              --->  coning at time n high_rate
% pos_ned              --->  Position Array (3x1) in LLH frame [rad] at time medium_rate -1
                                                                    
%Gravity (most time consuming part of ecef implementations)
Llh = pos_ned;

[Rn, Re, g, sL, cL, WIE_E]=geoparam(Llh);
Cne=pos2Cne(Llh(1), Llh(2));
ge=-Cne(:,3)*g;
                                                                          
                                                                        

%%Update attitude

%Part I:Nav Frame update
rot=-([0;0;WIE_E])*(dt);
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_b=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);
Cbe_ii =mx_b*Cbe;


%Part I:Body frame update
rot= rotation_angle+coning;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);
% Attitude update for ECEF frame rotation
Cbe_new = Cbe_ii*mx_a;


% Cbn_new=mx_b*Cbn*mx_a;

%%Update Velocity
%Body frame velocity rotation compensation
alfa_norm = norm( rotation_angle);
vr_a = (1/2)-(alfa_norm^2/24)+(alfa_norm^4/720);
vr_b = (1/6)-(alfa_norm^2/120)+(alfa_norm^4/5040);
v_rotation = vr_a*cross(rotation_angle,delta_vel)+vr_b*skew(rotation_angle)*cross(rotation_angle,delta_vel);
% Body integrated specific force increment
vel_b_specforce_increment = delta_vel+v_rotation+sculling;
%ECEF frame integrated specific force increment
vel_e_specforce_increment= Cbe_ii*vel_b_specforce_increment;
% integrated Coriolis acceleration and gravity increment
v_coriolis_gravity = (-ge+2*cross(Ve,[0;0;WIE_E]))*dt;
%ECEF frame velocity update
Ve_new = Ve + vel_e_specforce_increment+v_coriolis_gravity;

%% position update

ecef_new=ecef+0.5*(Ve+Ve_new)*dt;