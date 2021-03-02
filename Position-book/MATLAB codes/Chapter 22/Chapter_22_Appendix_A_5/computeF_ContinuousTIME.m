function phi = computeF_ContinuousTIME(r, F_e, Rb_e, tau_acc, tau_gyro, num_states)

%% INPUTS
% r            ---> last position array (3x1) computed with INS [m]
% F_e          ---> Accelerometer from body to Ecef frame matrix (3x3)
% Rb_e         ---> last Attitude matrix (3x3) computed from body to ECEF frame
% tau_acc      ---> Time constant array (3x1) of Gauss-Markov Accelerometers noise [hr]
% tau_gyro     ---> Time constant array (3x1) of Gauss-Markov Gyroscopes noise [hr]
% num-states   ---> Kalman Filter number of states

% f_e: accelerazioni nell'earth frame
% Re_b: matrice di rotazione da body ad earth
% alpha: parametro Gauss Markov giroscopi
% beta: parametro Gauss Markov accelerometri
% ordine: ordine di approssimazione di Taylor
global Ec;

F = zeros(num_states);

F(1:3,4:6)      = eye(3);
F(4:6,1:3)      = computeN(r); % gradiente gravitazionale
F(4:6,4:6)      = -2* skew([0;0;Ec.omega_e]);
F(4:6,7:9)      = -skew(F_e); 
F(4:6,10:12)    = Rb_e; 

F(7:9,7:9)      = -skew([0;0;Ec.omega_e]); 
F(7:9,13:15)    = Rb_e;

F(10:12,10:12)  = -diag(1./([tau_gyro(1), tau_gyro(2), tau_gyro(3)].*3600));
F(13:15,13:15)  = -diag(1./([tau_acc(1), tau_acc(2), tau_acc(3)].*3600)) ;
phi = F;
