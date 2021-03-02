function x = lls(X,r)
% LLS algorithm
% --------------------------------
% x = lls(X,r);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% 
L = size(X,2); % number of anchors
A = [-2*X' ones(L,1)];
b = r.^2-sum(X'.^2,2);
p = pinv(A'*A)*A'*b;
x= [p(1) ; p(2)];


