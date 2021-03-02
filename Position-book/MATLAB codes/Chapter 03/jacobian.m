function G = jacobian(X, x)
% Hessian matrix computation
% --------------------------------
% G = Jacobian matrix 
% G = jacobian(X, x)
% x = 2D position estimate
% X = Anchors position
%

[dim,L] = size(X); % L--number of anchors; dim--dimension of space
f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
G = (ones(L,1)*x' - X')./(f_TOA*ones(1,dim));
