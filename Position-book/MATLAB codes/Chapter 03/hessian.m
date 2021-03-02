function H = hessian(X,x,r)
% Hessian matrix computation
% --------------------------------
% H = hessiam(X,x,r);
% H = Hessian matrix 
% X = Anchors position
% x = 2D position estimate
%

L = size(X,2); % number of anchors
t1 = 0;
t2 = 0;
t3 =0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (x(1)-X(1,i))^2/ds(i)-(r(i)-ds(i)^(0.5))*(x(2)-X(2,i))^2/ds(i)^(1.5);
    t2 = t2 + (x(2)-X(2,i))^2/ds(i)-(r(i)-ds(i)^(0.5))*(x(1)-X(1,i))^2/ds(i)^(1.5);
    t3 = t3 + r(i)*(x(1)-X(1,i))*(x(2)-X(2,i))/ds(i)^(1.5);
end
H=2.*[t1 t3;
      t3 t2];
