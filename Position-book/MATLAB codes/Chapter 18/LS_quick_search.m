function y = LS_quick_search(p_cal,hat_bar_tau,P,M,B,c)

bar_tau_theory = (1/c) * sqrt( ( P(M+1:B,1) - p_cal(1,1) ).^2 + ( P(M+1:B,2) - p_cal(2,1) ).^2 );
y = sum((hat_bar_tau - bar_tau_theory).^2);

% bar_d_theory = sqrt( ( P(M+1:B,1) - p_cal(1,1) ).^2 + ( P(M+1:B,2) - p_cal(2,1) ).^2 );
% y = sum(( c * hat_bar_tau - bar_d_theory ).^2);