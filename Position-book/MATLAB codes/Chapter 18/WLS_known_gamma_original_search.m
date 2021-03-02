function y = WLS_known_gamma_original_search(p_cal,hat_bar_tau,bar_gamma,beta_bar,d_0,SNR,P,M,B,c,kappa)

bar_d_theory = sqrt( ( P(M+1:B,1) - p_cal(1,1) ).^2 + ( P(M+1:B,2) - p_cal(2,1) ).^2 );
bar_tau_theory = bar_d_theory / c;
fac = 1 + (bar_gamma.^2)./(16 * pi^2 * beta_bar^2 * bar_tau_theory.^2);
alpha_bar = sqrt( kappa * (d_0./(c*bar_tau_theory)).^(bar_gamma) ); 
sigma2_b = 1 ./ ( 8 * pi^2 * SNR * beta_bar^2 * alpha_bar.^2 .* fac);
y = sum( ((hat_bar_tau - bar_tau_theory).^2) ./ sigma2_b  );