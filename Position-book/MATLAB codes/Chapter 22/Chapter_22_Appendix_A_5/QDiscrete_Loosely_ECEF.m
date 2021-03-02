function parms_out = QDiscrete_Loosely_ECEF(dt_Execution, FDiscrete,G,white_noise_gyro, white_noise_acc,GM_gyro,GM_acc, tau_gyro, tau_acc)

%% INPUTS
% dt_Executions         ---> Time of KF execution [s]@ medium_rate 10Hz
% F_discrete            ---> F in discrete time matrix (15x15)
% G                     ---> G matrix in continuous time (5x4)
% white_noise_gyro      ---> Allan Variance array (3x1) of White Noise Gyros to convert as [rad/sqrt(sec)]^2
% white_noise_acc       ---> Allan Variance array (3x1) of White Noise Acc in [(m/s^2)/sqrt(sec)]^2
% GM_gyro               ---> Allan Variance array (3x1) of Gauss-Markov Noise Gyros in [rad/(sec)]
% GM_acc                ---> Allan Variance array (3x1) of Gauss-Markov Noise Acc in [(m/s^2)/(sec)]
% tau_gyro              ---> Allan Variance array (3x1) of Gauss-Markov time constant for Gyros [sec]
% tau_acc               ---> Allan Variance array (3x1) of Gauss-Markov time constant for Acc [sec]

dt             = dt_Execution;
fs             = 1/dt;

    [parms_in_Setting] = INS_Noise_parms(white_noise_gyro,white_noise_acc,GM_gyro,GM_acc,tau_gyro,tau_acc);
    Q = diag([...
        parms_in_Setting.q_a_white, parms_in_Setting.q_g_white,parms_in_Setting.q_b_a_GM,...
        parms_in_Setting.q_b_g_GM ]);
   
    
% Method to make Q a square matrix [ref 1 eq. 4.8]
    parms_out = 0.5 *(FDiscrete * G * Q * (G.') + ...
                  G * Q * (G.')* (FDiscrete.'))*dt;
%Qmatrix = parms_out;

end