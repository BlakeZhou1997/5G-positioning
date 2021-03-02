clear; clc;

Chapter_17_Example_3;

KVec = Rx_pos(:,1).^2 + Rx_pos(:,2).^2;
range_meas=recon_TOA*3e8;
hVec = (range_meas.^2-KVec);
GaMat = [-2*Rx_pos(:,1) -2*Rx_pos(:,2) ones(length(recon_TOA),1)];
theta0=[xHat yHat Rhat].';
[ZaVec_NP,ZaVec_NP_mat,fv_NP]=Chapter_17_Function_5(hVec,GaMat,theta0,[],5,1e-4);