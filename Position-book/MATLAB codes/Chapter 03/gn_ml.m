function x = gn_ml(X,r,iter,sigma2)
% Maximum Likelihood
% Local search scheme: Gauss-Netwon algorithm
% --------------------------------
% x = gn_ml(X,r,iter,sigma2);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
% sigma2 = variance of noise
% 
L = size(X,2); % number of anchors
x = lls(X,r);
for i = 1:iter
    G = jacobian(X, x);
    f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
    C_inv = diag(1./sigma2);
    x = x+pinv(G'*C_inv*G)*G'*C_inv*(r-f_TOA);
end