function H = basis_function(x)

r1 = x(1:3);
r2 = x(7:9);
b = norm(r2-r1);
doa12 = doa2(r2-r1);
psi = doa12(1);
phi = doa12(2);
dhtoadx = 1/b*[r2(1)-r1(1) r2(2)-r1(2) r2(3)-r1(3)];
dhdoadx = 1/b*[ -sin(psi)/cos(phi) cos(psi)/cos(phi) 0;
    -cos(psi)*sin(phi) -sin(psi)*sin(phi) cos(phi)];

r = norm(r1);
doar = doa2(r1);
lam = doar(1);
xi = doar(2);
dhradar = [r1(1)/r r1(2)/r r1(3)/r;
     -sin(lam)/(r*cos(xi)) cos(lam)/(r*cos(xi)) 0;
    -cos(lam)*sin(xi)/r -sin(lam)*sin(xi)/r cos(xi)/r];

H = [dhradar zeros(3,9);
    -dhtoadx zeros(1,3) dhtoadx zeros(1,3);
    -dhdoadx zeros(2,3) dhdoadx zeros(2,3)];
    



