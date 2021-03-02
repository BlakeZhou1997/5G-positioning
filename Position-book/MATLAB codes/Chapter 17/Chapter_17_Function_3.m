function [xHat,yHat,Rhat]=Chapter_17_Function_3(range_meas,Rx_pos)
% calculate the LS solution for geolocation based on
% range measurements
% INPUT:
%   range_meas
%   Rx_pos
% OUTPUT:
%   xHat
%   yHat

KVec = Rx_pos(:,1).^2 + Rx_pos(:,2).^2;
hVec = (range_meas.^2-KVec);
GaMat = [-2*Rx_pos(:,1) -2*Rx_pos(:,2) ones(length(range_meas),1)];
ZaVec = inv(GaMat'*GaMat)*GaMat'*hVec;
xHat = ZaVec(1); yHat = ZaVec(2); Rhat = ZaVec(3);