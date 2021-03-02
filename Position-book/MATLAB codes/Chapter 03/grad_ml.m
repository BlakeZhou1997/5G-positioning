function g = grad_ml(X,x,r,sigma2)
% Gradient computation in Maximun Likelihood
% --------------------------------
% g = grad(X,x,r,sigma2);
% g = gradient vector 
% X = Anchors position
% x = 2D position estimate
% r = TOA measurement vector
% sigma2 = variance of noise
%
L = size(X,2); % number of anchors
t1 = 0;
t2 = 0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (1/sigma2(i))*(r(i)-ds(i)^(0.5))*(x(1)-X(1,i))/ds(i)^(0.5);
    t2 = t2 + (1/sigma2(i))*(r(i)-ds(i)^(0.5))*(x(2)-X(2,i))/ds(i)^(0.5);
end
g=-2.*[t1; t2];