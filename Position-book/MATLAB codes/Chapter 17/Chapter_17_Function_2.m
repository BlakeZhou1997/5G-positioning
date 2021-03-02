function [trim_indices] = ...
    Chapter_17_Function_2(meas_toa,Porder,alpha,tail);
% Trimmed-mean algorithm with a polynomial fit to ToA measurements
% INPUT:
%   meas_toa    : toa measurements vector
%   Porder      : polynomial order
%   alpha       : parameter for Gaussianity SW test
%   tail        : parameter for Gaussianity SW test
% OUTPUT:
%   trim_indices: indices of trimmed TOA

Ktheta=length(meas_toa);
poly_coeff = polyfit(1:Ktheta,meas_toa.',Porder); 
poly_fit = polyval(poly_coeff,1:Ktheta);
residualTOA = meas_toa - poly_fit.'; 
[sort_TOA,trim_idx_TOA] = sort(residualTOA,'ascend');

%% trimming iteration 
H_TOA = 1; % flag for iteration
counter_trim_TOA = 0; % init
trimmer_per_TOA = 5; % no of measurements to be trimmed every iteration
trim_per_init = 5;
while (H_TOA)
    trim_TOA = sort_TOA(1:Ktheta-trimmer_per_TOA); % sorting
    [H_TOA, pValue, W] = Chapter_17_swtest(trim_TOA,alpha,tail); % SWtest
    counter_trim_TOA = counter_trim_TOA + 1; % increase counter

    if counter_trim_TOA >= 5 % iteration limit
        H_TOA = 0; break
    end
    if H_TOA == 1
        trimmer_per_TOA = trim_per_init + ...
            counter_trim_TOA*floor(Ktheta/10); 
    end
end
trim_indices = trim_idx_TOA(Ktheta-trimmer_per_TOA+1:Ktheta);
