% RMSE of the position estimate as a function of the transmitted SNR
clc
clear all
c = 3 * 10^8;
beta_bar = 2 * pi * (1/sqrt(3)) * 5 * 10^6;
r = 2000; %4000/sqrt(3);
p = (1/2) * r * cos(pi/6)* [cos(pi/6);sin(pi/6)];
B = 7;
M = 0; % base stations receiving NLOS
a_range = [4.6 4.0 3.6];
b_range = [0.0075 0.0065 0.0050];
c_range = [12.6 17.1 20.0];
scenario_range = [1:length(a_range)];
p_true = p;
h_b = 45; % (80+10)/2c
f = 1.9 * 10^9;
d_0 = 100;
kappa = c^2 / (16 * pi^2 * f^2 * d_0^2);
P = r * [0,3/2,0,-3/2,-3/2,0,3/2;0,sqrt(3)/2,sqrt(3),sqrt(3)/2,-sqrt(3)/2,-sqrt(3),-sqrt(3)/2].';
P_diff = P - repmat(p',B,1);
phi = atan2(P_diff(:,2),P_diff(:,1));
Phi = [cos(phi).';sin(phi).'];
Phi_tilde = Phi(:,1:M);
Phi_bar = Phi(:,M+1:B);
n_scenario = 0;
N_R = 10000;
SNR_range = [50:10:230]; %dB
n_SNR = 0;
for SNR_dB = SNR_range;
    n_SNR = n_SNR + 1;
    SNR = 10^(SNR_dB/10);
    for scenario = 1  %scenario_range
        n_scenario = n_scenario + 1;
        gamma_gen = a_range(scenario) - b_range(scenario)*h_b + c_range(scenario)/h_b;
        gamma_true = ones(B,1)*gamma_gen;
        bar_gamma = gamma_true(M+1:B,1);
        for b = 1 : B
                x_tilde(b) = P(b,1) - p_true(1);
                y_tilde(b) = P(b,2) - p_true(2);
                if b <= M    
                    tau_true(b) = (1/c) * (sqrt(x_tilde(b)^2 + y_tilde(b)^2) + l_true(b));
                else
                    tau_true(b) = (1/c) * (sqrt(x_tilde(b)^2 + y_tilde(b)^2));
                end
                d_true(b) = c * tau_true(b);
                alpha_true(b) = sqrt( kappa * ( d_0/d_true(b) )^(gamma_true(b)) );
                frac(b) = 1 + (gamma_true(b)^2)/(8 * pi^2 * beta_bar^2 * tau_true(b)^2);
                sigma2_true(b) = 1 / ( 8 * pi^2 * SNR * alpha_true(b)^2 * beta_bar^2 * frac(b));
        end % b
        bar_sigma2 = sigma2_true(M+1:B)';
        for n_R = 1 : N_R
            for b = 1 : B
                u_1(b) = randn;
                hat_tau_gen(b) = tau_true(b) + sqrt(sigma2_true(b)) * u_1(b);
            end
            hat_bar_tau = hat_tau_gen(M+1:B)';
            p_in = p_true + randn(1,1);
            bar_eta_in = [p_in;gamma_true(M+1,1)] + randn(3,1);
            hat_p_LS = fminsearch('LS_quick_search',p_in,[],hat_bar_tau,P,M,B,c);
            hat_p_WLS_known_gamma = fminsearch('WLS_known_gamma_original_search',p_in,[],hat_bar_tau,bar_gamma,beta_bar,d_0,SNR,P,M,B,c,kappa);
            hat_p_ML_known_gamma  =  fminsearch('ML_known_gamma_original_search',p_in,[],hat_bar_tau,bar_gamma,beta_bar,d_0,SNR,P,M,B,c,kappa);
            SE_p_LS(n_SNR,n_R) = sum((hat_p_LS - p_true).^2);
            SE_p_WLS_known_gamma(n_SNR,n_R) = sum((hat_p_WLS_known_gamma - p_true).^2);
            SE_p_ML_known_gamma(n_SNR,n_R) = sum((hat_p_ML_known_gamma - p_true).^2);
            error_p_LS(:,n_R) = hat_p_LS - p_true;
            error_p_WLS_known_gamma(:,n_R) = hat_p_WLS_known_gamma - p_true;
            error_p_ML_known_gamma(:,n_R) = hat_p_ML_known_gamma - p_true;
        end
        RMSE_p_LS(n_SNR) = sqrt(mean(SE_p_LS(n_SNR,:)))
        RMSE_p_WLS_known_gamma(n_SNR) = sqrt(mean(SE_p_WLS_known_gamma(n_SNR,:)))
        B_LS_theory = c^2 * inv(Phi_bar * Phi_bar') * Phi_bar * diag( sigma2_true(M+1:B) ) * Phi_bar' * inv(Phi_bar * Phi_bar');
        RMSE_p_LS_theory(n_SNR) = sqrt( trace( B_LS_theory ) );
        RMSE_p_ML_known_gamma(n_SNR) = sqrt(mean(SE_p_ML_known_gamma(n_SNR,:)))
        B_pp_given_gamma =  c^2 * (1 / ( 8 * pi^2 * beta_bar^2 * SNR)) *  inv(Phi_bar* diag(alpha_true(M+1:B))^2 * (eye(B-M) + ...
                            (1/(16 * pi^2 * beta_bar^2)) * diag(gamma_true(M+1:B))^2 * inv( diag(tau_true(M+1:B))^2 ) ) * Phi_bar');
        bias_p_LS(n_SNR) = mean(sqrt( error_p_LS(1,:).^2 + error_p_LS(2,:).^2 ));
        bias_p_WLS_known_gamma(n_SNR) =  mean(sqrt(error_p_WLS_known_gamma(1,:).^2 + error_p_WLS_known_gamma(2,:).^2 ));
        bias_p_ML_known_gamma(n_SNR) = mean(sqrt( error_p_ML_known_gamma(1,:).^2 + error_p_ML_known_gamma(2,:).^2 ));
        CRB_p_given_gamma(n_SNR) = sqrt(trace(B_pp_given_gamma));
    end
end

 figure
 hold on
 %subplot(2,1,1);
 grid on
 plot(SNR_range,RMSE_p_LS,'r*','LineWidth',.5,'MarkerEdgeColor','r','MarkerSize',5)
 plot(SNR_range,RMSE_p_LS_theory,':r','LineWidth',.5,'MarkerEdgeColor','r','MarkerSize',5)
 plot(SNR_range,RMSE_p_WLS_known_gamma,'mo','LineWidth',.5,'MarkerEdgeColor','m','MarkerSize',5)
 plot(SNR_range,RMSE_p_ML_known_gamma,'bs','LineWidth',.5,'MarkerEdgeColor','b','MarkerSize',5)
 plot(SNR_range,CRB_p_given_gamma,'-k','LineWidth',.5,'MarkerEdgeColor','k','MarkerSize',5)
 legend('LS: simulation','LS: analysis','WLS: simulation','ML: simulation','CRB')
 h = findobj(gcf,'type','axes','tag','legend');
 Pos = get(h,'position');
 Pos(3) = 1.9 * Pos(3); % Double the length
 Pos(4) = 3 * Pos(4); % Double the length
 set(h,'position',Pos) % Implement it
 title('RMSE as a function of SNR ($M=0$, $N_R=1{,}000$)')
 ylabel('Root Mean Square Error (m)')%,'FontSize',10,'FontName','Times New Roman')
 xlabel('Transmitted SNR $\frac{E_{\mathrm{s}}}{\sigma^2_{\mathrm{n}}}$ (dB)')%,'FontSize',10,'FontName','Times New Roman')
 %axis([min(xi_range) max(xi_range) min(varepsilon_CRB(4,:)) max(varepsilon_CRB(3,:))])
 hold off
 
 figure
 hold on
 %subplot(2,1,2);
 grid on
 plot(SNR_range,bias_p_LS,'r*','LineWidth',.5,'MarkerEdgeColor','r','MarkerSize',5)
 plot(SNR_range,bias_p_WLS_known_gamma,'mo','LineWidth',.5,'MarkerEdgeColor','m','MarkerSize',5)
 plot(SNR_range,bias_p_ML_known_gamma,'bs','LineWidth',.5,'MarkerEdgeColor','b','MarkerSize',5)
 legend('LS: simulation','WLS: simulation','ML: simulation','CRB')
 h = findobj(gcf,'type','axes','tag','legend');
 Pos = get(h,'position');
 Pos(3) = 1.9 * Pos(3); % Double the length
 Pos(4) = 3 * Pos(4); % Double the length
 set(h,'position',Pos) % Implement it
 title('Bias as a function of SNR ($M=0$, $N_R=1{,}000$)')
 ylabel('Bias (m)')%,'FontSize',10,'FontName','Times New Roman')
 xlabel('Transmitted SNR $\frac{E_{\mathrm{s}}}{\sigma^2_{\mathrm{n}}}$ (dB)')%,'FontSize',10,'FontName','Times New Roman')
 %axis([min(xi_range) max(xi_range) min(varepsilon_CRB(4,:)) max(varepsilon_CRB(3,:))])
 hold off
 
 set(0,'defaulttextinterpreter','none')%laprint