function parms_out = FDiscrete_ECEF(dt_INS,F)

dt             = dt_INS;

parms_out = eye(length(F(1,:))) + F*dt + ((F*dt)^2)/2;
 
end
