function x = nr_ml(X,r,iter,sigma2)
% Maximum Likelihood
% Local search scheme: Newton-Raphson algorithm
% --------------------------------
% x = nr_ml(X,r,iter,sigma2);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
% sigma2 = variance of noise
%

x=lls(X,r);
for i=1:iter
    H=hessian_ml(X,x,r,sigma2);
    g=grad_ml(X,x,r,sigma2);
    x=x-pinv(H)*g;

end



