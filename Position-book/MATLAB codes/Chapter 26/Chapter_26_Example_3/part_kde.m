% Compute KDE (y) at given points x (i.e., compute weights of the particles or grid-points)
% Gaussian spherical kernel used
function y=part_kde(wei,mean,band,x);
% (wei, mean) - weighted particles, band - bandwidth, x - desired points
epsilon=0.000001; % some very small number
y=zeros(size(x));
for k=1:length(wei)
    % log of one component
    log_exp=0.5*abs( (real(x-mean(k)))/real(band)+j*(imag(x-mean(k)))/imag(band) ).^2;
    min_log_exp=min(min(log_exp));
    y=y+wei(k)*exp(-(log_exp-min_log_exp));
end
y=y/sum(sum(y)); % normalize (y can be 2D matrix or vector)
if round(sum(sum(y)))==1 % always, except if equal to NaN
    y=y+epsilon; % avoid zero
    y=y/sum(sum(y));  % normalize again
else % never should happen
    y=ones(size(x)); y=y/sum(sum(y)); 
    fprintf('Warning: undefined result (NaN) overrided using uniform distribution\n');
end  
end 

