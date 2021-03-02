function x = nr_nls(X,r,iter)
% Nonlinear Least Squares
% Local search scheme: Newton-Raphson algorithm
% --------------------------------
% x = nr_nls(X,r,iter);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
%

x=lls(X,r);
for i=1:iter
    H=hessian(X,x,r);
    g=grad(X,x,r);
    x=x-pinv(H)*g;
end



