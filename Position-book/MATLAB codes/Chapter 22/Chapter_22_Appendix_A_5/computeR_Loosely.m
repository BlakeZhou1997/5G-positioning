function [R]=computeR_Loosely()

posSigma2 = 1e1; % m^2, variance of GPS position X-Y-Z [m^2]
velSigma2 = 1e-1; % m^2/s^2, variance of GPS velocity X-Y-Z {m/s]^2

vect_pos = [posSigma2, posSigma2,posSigma2];
vect_vel = [velSigma2, velSigma2,velSigma2];

R=diag([vect_pos, vect_vel]);

return;