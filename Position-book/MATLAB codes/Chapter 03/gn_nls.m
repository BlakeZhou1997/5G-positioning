function x = gn_nls(X,r,iter)
% Nonlinear Least Squares
% Local search scheme: Gauss-Netwon algorithm
% --------------------------------
% x = gn_nls(X,r,iter);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% iter = number of iteration
%
L = size(X,2); % number of anchors
x = lls(X,r);
for i = 1:iter
    G = jacobian(X, x);
    f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
    x = x+pinv(G'*G)*G'*(r-f_TOA);
end