function KalmanStates = LooselyCoupling_ECEF_ASSIST(INS_pos_ECEF,GPS_pos_ECEF, INS_vel_ECEF, GPS_vel_ECEF, Num_states,P,QDisc)

%% INPUTS

% INS_pos_ECEF      ---> INS position array (3x1) in ECEF
% GPS_pos_ECEF      ---> GPS position array (3x1) in ECEF
% INS_vel_ECEF      ---> INS velocity array (3x1) in ECEF
% GPS_vel_ECEF      ---> GPS velocity array (3x1) in ECEF
% Num_states        ---> Kalman filter number of states
% P                 ---> Kalman filter P error matrix (15x15)
% QDisc             ---> Discrete Matrix of Kalman Model 


%% Computing the INS-GPS misclosure
    
zCurrentMeasurement = [GPS_pos_ECEF; GPS_vel_ECEF];
NominalMeasurement = [INS_pos_ECEF;INS_vel_ECEF];

KFMeasurement = zCurrentMeasurement - NominalMeasurement ;
                                                                    

%% computing the Transition Matrix

HMatrix =  [  1                      zeros(1,Num_states-1);
              0   1                  zeros(1,Num_states-2);
              0   0   1              zeros(1,Num_states-3);
              0   0   0   1          zeros(1,Num_states-4);       
              0   0   0   0  1       zeros(1,Num_states-5);
              0   0   0   0  0  1    zeros(1,Num_states-6)];
%% computing the noise Measurement covariance Matrix

[R]=computeR_Loosely();

%% Computing the Kalman filter 
xn_1 = zeros(Num_states,1); % zero previous (error) state 
KalmanOutputs = ComputeKalman2(xn_1,P,QDisc,HMatrix,R,KFMeasurement); % JTM 3/3/11
KalmanStates.P = KalmanOutputs.P;
KalmanStates.X = KalmanOutputs.X;
end