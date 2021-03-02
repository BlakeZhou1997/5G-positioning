
function parms_out = INS_Noise_parms(white_noise_gyro,white_noise_acc,GM_gyro,GM_acc,tau_gyro,tau_acc)

%% INPUT
% white_noise_gyro      ---> Allan Variance array (3x1) of White Noise Gyros in [deg/sqrt(hr)]
% white_noise_acc       ---> Allan Variance array (3x1) of White Noise Acc in [G/sqrt(hr)]
% GM_gyro               ---> Allan Variance array (3x1) of Gauss-Markov Noise Gyros in [deg/(hr)]
% GM_acc                ---> Allan Variance array (3x1) of Gauss-Markov Noise Acc in [G/(hr)]
% tau_gyro              ---> Allan Variance array (3x1) of Gauss-Markov time constant for Gyros [hr]
% tau_acc               ---> Allan Variance array (3x1) of Gauss-Markov time constant for Acc [hr]

%%computed by Allan Variance Technique and Reference D2.2


%% Real values mine
%ARW 
sdgyr_noise =  [white_noise_gyro(1), white_noise_gyro(2), white_noise_gyro(3)]*pi/(180*sqrt(3600)); 
sdacc_noise = [white_noise_acc(1),white_noise_acc(2),white_noise_acc(3)].*(9.80655/sqrt(3600));

parms_out.q_g_white = 1e0.*sdgyr_noise.^2 ;
parms_out.q_a_white = 1e0.*sdacc_noise.^2;

%Bias Instability-Flicker Noise-GM % see Yuksel document for reference
Tba_1 = [tau_acc(1)/1.89, tau_acc(2)/1.89, tau_acc(3)/1.89].*3600; 
sd_ba1 = [GM_acc(1), GM_acc(2) , GM_acc(3)].*9.80655;  
sd_ba1 = sd_ba1./(0.437*sqrt(Tba_1));                  %qc noise amplitude 
% parms_out.q_b_a_GM =1e0.*sd_ba1.^2;
parms_out.q_b_a_GM =1e0.*[1e0  1e0  1e3].*sd_ba1;


Tbg_1 = [tau_gyro(1)/1.89, tau_gyro(2)/1.89, tau_gyro(3)/1.89].*3600;
sd_bg1 =  [GM_gyro(1), GM_gyro(2) , GM_gyro(3)].*(pi/180)./3600; 
sd_bg1 = sd_bg1./(0.437*sqrt(Tbg_1));                  %qc noise amplitude 
parms_out.q_b_g_GM =1e0.*sd_bg1;
%  parms_out.q_b_g_GM =1e0.*[1e0 5e2 5e-1].*sd_bg1.^2;

end