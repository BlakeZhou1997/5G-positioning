function x = sd_nls(X,r,iter,mu)
% Nonlinear Least Squares
% Local search scheme: Steepest descent algorithm
% --------------------------------
% x = sd_nls(X,r,iter,mu);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
% mu = step size
%

x=lls(X,r);
for i=1:iter
    g=grad(X,x,r);
    x=x-mu*g;
end

