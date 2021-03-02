function pdot = covprop(t,ptemp,f,g,q)

lp = length(ptemp);
p = reshape(ptemp,sqrt(lp),sqrt(lp));

pdot = f*p+p*f'+g*q*g';

ptemp = pdot;
pdot = reshape(ptemp,lp,1);