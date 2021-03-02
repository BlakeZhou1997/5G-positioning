clear; clc;

Chapter_17_Example_1;

Porder=3;
alpha=0.05;
tail=0;
[trim_indices]=Chapter_17_Function_2(meas_toa,Porder,alpha,tail);

trimmed_TOA = meas_toa;
trimmed_TOA(trim_indices) = 0;
Rx_pos(trim_indices,:)=0;
[xHat,yHat,Rhat]=Chapter_17_Function_3(trimmed_TOA*3e8,Rx_pos);