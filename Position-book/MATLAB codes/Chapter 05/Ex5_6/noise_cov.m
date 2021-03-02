function R = noise_cov(x,sigr,sigt,sigp)

r = x(1);
p = x(2);
t = x(3);


Rk = diag([(r^2+sigr^2)*(1+exp(-2*sigt^2))*(1+exp(-2*sigp^2))/4-r^2*exp(-(sigt^2+sigp^2))
    (r^2+sigr^2)*(1-exp(-2*sigt^2))*(1+exp(-2*sigp^2))/4
    (r^2+sigr^2)*(1+exp(-2*sigt^2))*(1-exp(-2*sigp^2))/4
    (r^2+sigr^2)*(1-exp(-2*sigt^2))/2]);

Ck = [cos(p)*cos(t) -sin(t)*cos(p) -cos(t)*sin(p);
    sin(t)*cos(p) cos(t)*cos(p) -sin(t)*sin(p);
    sin(p) 0 cos(p)];

Tk = [Ck [sin(t)*sin(p); -cos(t)*sin(p);0]];

R = Tk*Rk*Tk';