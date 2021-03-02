% TOA Based Positioning
% main function

close all
clear all
clc

x1 = [0,0];
x2 = [0,10];
x3 = [10,0];
x4 = [10,10];
X = [x1;x2;x3;x4]'; % anchors position
L = size(X,2); % number of anchors
x = [2,3]';
d = sqrt(sum((x*ones(1,L)-X).^2,1)); 
d = d';
not = 1e2; % number of trails
range_dB =-10:5:60; % range of SNR (in dB)
m = 0;
for dB = range_dB
    sigma2 = d.^2/10^(dB/10); % sigma2--square of sigma, here we use: SNR_dB = 10log(d^2/sigma^2)
    for run = 1:not
        r = d + randn(L,1).*sqrt(sigma2);
        % Linear Approaches
        x_lls(run,:) = lls(X,r); 
        x_wls(run,:) = wls(X,r,sigma2);
        x_wls2(run,:) = wls2(X,r,sigma2);
        x_sub(run,:) = sub(X,r);
        
        % Nonlinear Approaches
        iter = 100;
        % == Newton-Raphson algorithm == %
        x_nr_nls(run,:) = nr_nls(X,r,iter);
        x_nr_ml(run,:) = nr_ml(X,r,iter,sigma2);
        % == Gauss-Netwon algorithm == %
        x_gn_nls(run,:) = gn_nls(X,r,iter);
        x_gn_ml(run,:) = gn_ml(X,r,iter,sigma2);
        % == Steepest descent algorithm == %
        mu_nls = 0.1;
        x_sd_nls(run,:) = sd_nls(X,r,iter,mu_nls);
        mu_ml = 0.00001;
        x_sd_ml(run,:) = sd_ml(X,r,iter,mu_ml,sigma2);
    end;
    m = m + 1
    mse_lls(m) = mean(sum((x_lls - ones(not, 1)*x').^2, 2));
    mse_wls(m) = mean(sum((x_wls - ones(not, 1)*x').^2, 2));
    mse_wls2(m) = mean(sum((x_wls2 - ones(not, 1)*x').^2, 2));
    mse_sub(m) = mean(sum((x_sub - ones(not, 1)*x').^2, 2));
    
    mse_nr_nls(m) = mean(sum((x_nr_nls - ones(not, 1)*x').^2, 2));
    mse_nr_ml(m) = mean(sum((x_nr_ml - ones(not, 1)*x').^2, 2));
    
    mse_gn_nls(m) = mean(sum((x_gn_nls - ones(not, 1)*x').^2, 2));
    mse_gn_ml(m) = mean(sum((x_gn_ml - ones(not, 1)*x').^2, 2));
    
    mse_sd_nls(m) = mean(sum((x_sd_nls - ones(not, 1)*x').^2, 2));
    mse_sd_ml(m) = mean(sum((x_sd_ml - ones(not, 1)*x').^2, 2));

    crlb(m) = CRLB([x X]', sigma2);
end

figure
plot(X(1,:), X(2,:), 'bs'), grid on, hold on
plot(x(1),x(2), 'ro'), grid on, hold on
legend('Anchors', 'Unknown-position sensor')
xlabel('x-coordinate (m)')
ylabel('y-coordinate (m)')


figure
plot(range_dB, 10*log10(mse_lls), 'r.', range_dB,10*log10(mse_wls), 'bo',...
    range_dB,10*log10(mse_wls2), 'm*',range_dB, 10*log10(mse_sub), 'cd',...
    range_dB, 10*log10(crlb), 'k-');
legend('LLS','WLS','Two-step WLS','Subspace','CRLB');

figure;
plot(range_dB, 10*log10(mse_nr_nls), 'bo', range_dB,10*log10(mse_nr_ml), 'r*',range_dB, 10*log10(crlb), 'k-');
legend('NR-NLS','NR-ML','CRLB');

% figure;
% plot(range_dB, 10*log10(mse_gn_nls), 'bo', range_dB,10*log10(mse_gn_ml), 'r*',range_dB, 10*log10(crlb), 'k-');
% legend('GN-NLS','GN-ML','CRLB');
% 
% figure;
% plot(range_dB, 10*log10(mse_sd_nls), 'bo', range_dB,10*log10(mse_sd_ml), 'r*',range_dB, 10*log10(crlb), 'k-');
% legend('SD-NLS','SD-ML','CRLB');

