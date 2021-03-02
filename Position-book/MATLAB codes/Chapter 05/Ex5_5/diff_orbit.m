function h=diff_orbit(R,mu)

r1=R(1);
r2=R(2);
r3=R(3);
h=zeros(3);


S=mu/norm(R)^5;

h=S*[(2*r1^2-r2^2-r3^2) 3*r1*r2 3*r1*r3;
        3*r1*r2 (2*r2^2-r1^2-r3^2) 3*r2*r3;
        3*r1*r3 3*r2*r3 (2*r3^2-r1^2-r2^2);];
