function x = wls2(X,r,sigma2)
% Two-step Weighted LLS algorithm
% --------------------------------
% x = wls2(X,r,sigma2);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% sigma2 = variance of noise
% 
L = size(X,2); % number of anchors

% The first step
A = [-2*X' ones(L,1)];
b = r.^2-sum(X'.^2,2);
W = 1/4*diag(1./(sigma2.*r.^2));
z = pinv(A'*W*A)*A'*W*b;

% The second step
s = sign(z(1:2));
G = [1 0;0 1;1 1];
h = [z(1)^2;z(2)^2;z(3)];
Phi = diag(2*[z(1:2);1])*pinv(A'*W*A)*diag(2*[z(1:2);1]);
z = pinv(G'*Phi*G)*G'*Phi*h;
x = real(sign(s).*sqrt(z));

