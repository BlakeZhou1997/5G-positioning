function T = find_delay(x,c)

r12 = x(1:3)';
v2 = x(4:6)';

A = c^2-v2'*v2;
B = 2*v2'*r12;
C = -r12'*r12;

T = (-B+sqrt(B^2-4*A*C))/2/A;