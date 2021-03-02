function y=orbit(t,x,mu)
y = zeros(6,1);
y(1:3)=x(4:6);
y(4:6) = -mu/norm(x(1:3))^3*x(1:3);
%y=y';