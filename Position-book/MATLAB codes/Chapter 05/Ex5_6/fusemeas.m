function [yf,Rf,h,Rcm,ycm]=fusemeas(ym,tof,Rm,sr,sp,c,tk)
% Fuse all TOA and DOA measurements based on TOF and noise variance
% Required matlab files
%      noise_cov.m, polar2car.m

Rcm = [];
ycm = [];

len = length(ym);
y = reshape(ym,3,len/3);
y(1,:) = (y(1,:))/2*c; % divide the TOF by 2 = approximate range between UAV and beacon


sum_w = zeros(3);
yf = zeros(3,1);
[row,~]=size(Rm);
Re = kron(Rm(tk,:),ones(row,1))-Rm;
h = zeros(3,6);


% Calculate the weighting factor based on TOF
if length(tof) == 1  %add to avoid singularity when determine the weighting factor
    weight_vec = 1;
    weight_rvec = 1;
else
    DT = tof(tk)-tof;
    DT = 1 - DT/norm(DT);
    weight_vec = DT.^2/(DT'*DT);
end

for jj=1:len/3
    Rt = noise_cov(y(:,jj),sr,sp,sp); %Compute covariance in cartersian
    yc = polar2car(y(:,jj)); % convert TOA, DOA into cartesian
    Rcm = [Rcm;Rt];    
    ycm = [ycm;yc];
    
    % The summation of inverse covariance weight fusion
    Rtinv = inv(Rt);
    yf = yf+weight_vec(jj)*Rtinv*(Re(jj,:)'+yc);   
    sum_w = sum_w+weight_vec(jj)*Rtinv;    
    h = h+weight_vec(jj)*Rtinv*([-eye(3) zeros(3)]);
end

% The multiplication of covariance and summation of inverse covariance
% weighting
Rf = inv(sum_w);
yf = Rf*yf; 
h = Rf*h;
Rf = length(tof)*Rf;