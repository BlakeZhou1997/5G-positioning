function [recon_TOA] = Chapter_17_Function_4(meas_toa,...
    Rx_pos,Porder,alpha,tail)
% Iterative algorithm that re-construct the TOA measurement
% based on geometrical estimation from LS solution
% INPUT:
%   meas_toa        : TOA measurements
%   Rx_pos          : Positions where the TOA measurements 
%                     are collected
%   Porder          : polynomial order
%   (alpha,tail)    : parameters for Gaussianity SW test
% OUTPUT:
%   recon_TOA       : re-constructed TOA measurements


Rx_x1 = Rx_pos(:,1);
Rx_y1 = Rx_pos(:,2);
c=3e8; 

%% Step 1: Trim TOA measurements
[trim_indices] = Chapter_17_Function_2(meas_toa,Porder,alpha,tail);

%% Step 2: Estimate MT position from the remaining TOA
trimmed_TOA = meas_toa; trimmed_TOA(trim_indices) = 0;
rVec_trim = trimmed_TOA*3e8; % range measurement
KVec_trim = Rx_x1.^2 + Rx_y1.^2; 
KVec_trim(trim_indices) = 0;
hVec_trim = (rVec_trim.^2-KVec_trim); 
GaMat_trim = [-2*Rx_x1 -2*Rx_y1 ones(length(meas_toa),1)];
GaMat_trim(trim_indices,:) = 0;
ZaVec_trim = inv(GaMat_trim'*GaMat_trim)*GaMat_trim'...
    *hVec_trim;

%% Step 3: Iteration
% initialize
trimmed_TOA01 = trimmed_TOA; 
trim_indices01 = trim_indices;
estTxX_cur = ZaVec_trim(1); 
estTxY_cur = ZaVec_trim(2);
tolerance_min = 100; tolerance_max = 1000; 
ZaVec_trim_old = ZaVec_trim;

loop_cnt = 1; count = 1; 
while (loop_cnt)
    
    trimming_toa = trimmed_TOA01;
    % replace the trimmed TOA based on prev estimated 
    % MT position
    for kk = sort(trim_indices01)
        trimming_toa(kk) = 1/c*sqrt((ZaVec_trim_old(1)-...
            Rx_x1(kk)).^2 + (ZaVec_trim_old(2)-...
            Rx_y1(kk)).^2);
    end
    
    % trim iteration
    [trim_indices01] = Chapter_17_Function_2(trimming_toa,Porder,alpha,tail);
    
    % estimate MT position from the remaining TOA
    trimmed_TOA01 = trimming_toa; 
    trimmed_TOA01(trim_indices01) = 0;
    rVec_trim01 = trimmed_TOA01*3e8; % range measurement
    KVec_trim01 = Rx_x1.^2 + Rx_y1.^2; 
    KVec_trim01(trim_indices01) = 0;
    hVec_trim01 = (rVec_trim01.^2-KVec_trim01);
    GaMat_trim01 = [-2*Rx_x1 -2*Rx_y1 ones(length(meas_toa),1)];
    GaMat_trim(trim_indices01,:) = 0;
    ZaVec_trim01 = inv(GaMat_trim01'*GaMat_trim01)...
        *GaMat_trim01'*hVec_trim01;
    
    % check stopping criterion
    err_new = sqrt((estTxX_cur - ZaVec_trim01(1)).^2 ...
        + (estTxY_cur - ZaVec_trim01(2)).^2);
    if err_new <= tolerance_min 
        ZaVec_trimming = ZaVec_trim01;
        loop_cnt = 0; break;
    elseif err_new >= tolerance_max
        ZaVec_trimming = ZaVec_trim_old;
        loop_cnt = 0; break;
    elseif count >= 5
        ZaVec_trimming = ZaVec_trim01;
        loop_cnt = 0; break;
    else
        count = count + 1;
        estTxX_cur = ZaVec_trim01(1); 
        estTxY_cur = ZaVec_trim01(2);
        ZaVec_trim_old = ZaVec_trim01;
    end
end

%% Step 4: Reconstruct TOA based on last iteration
recon_TOA = trimmed_TOA01;
% replace the trimmed TOA based on prev estimated 
% MT position
for kk = sort(trim_indices01)
    recon_TOA(kk) = 1/c*sqrt((ZaVec_trimming(1)-...
        Rx_x1(kk)).^2 + (ZaVec_trimming(2)-Rx_y1(kk)).^2);
end