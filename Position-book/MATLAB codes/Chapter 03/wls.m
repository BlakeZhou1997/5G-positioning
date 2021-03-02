function x = wls(X,r,sigma2)
% Weighted LLS algorithm
% --------------------------------
% x = wls(X,r,sigma2);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% sigma2 = variance of noise
% 

L = size(X,2); % number of anchors
A = [-2*X' ones(L,1)];
b = r.^2-sum(X'.^2,2);
W = 1/4*diag(1./(sigma2.*r.^2));
p = pinv(A'*W*A)*A'*W*b;
x= [p(1) ; p(2)];


