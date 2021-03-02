function y = polar2car(x)

y=zeros(3,1);
y(1) = x(1)*cos(x(2))*cos(x(3));
y(2) = x(1)*sin(x(2))*cos(x(3));
y(3) = x(1)*sin(x(3));
