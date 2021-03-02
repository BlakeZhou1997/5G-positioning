function x = sd_ml(X,r,iter,mu,sigma2)
% Maximum Likelihood
% Local search scheme: Steepest descent algorithm
% --------------------------------
% x = sd_nls(X,r,iter,mu,sigma2);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
% mu = step size
% sigma2 = variance of noise
%

x=lls(X,r);
for i=1:iter
    g=grad_ml(X,x,r,sigma2);
    x=x-mu*g;
end
