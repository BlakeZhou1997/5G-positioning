function x = sub(X,r)
% Subspace algorithm
% --------------------------------
% x = sub(X,r);
% x = 2D position estimate
% X = Anchors position
% r = TOA measurement vector
% 
Y = X';
L = size(Y,1); % number of anchors
R = squareform(pdist(Y));
D = zeros(L);
for i=1:L
    for j=1:L
        D(i,j)=0.5*(r(i)^2+r(j)^2-R(i,j)^2);
    end
end
% [U,Lamda] = eig(D);
[U,S,V] = svd(D);
Un = U(:,3:end);
x = (Y'*(Un*Un')*ones(L,1))/(ones(L,1)'*(Un*Un')*ones(L,1));
x = x';