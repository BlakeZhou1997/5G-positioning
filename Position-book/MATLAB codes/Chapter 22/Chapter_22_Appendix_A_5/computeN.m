function N = computeN(r)



global Const;

R = sqrt(sum(r.^2));

N = zeros(3);

% N(1,1) = (Ec.mu/(R^3))*(3*r(1)/R^2-1)+Ec.omega_e^2;
% N(1,2) = (Ec.mu/(R^3))*(3*r(1)*r(2)/R^2);
% N(1,3) = (Ec.mu/(R^3))*(3*r(1)*r(3)/R^2);
% N(2,1) = (Ec.mu/(R^3))*(3*r(1)*r(2)/R^2);
% N(2,2) = (Ec.mu/(R^3))*(3*r(2)/R^2-1)+Ec.omega_e^2;
% N(2,3) = (Ec.mu/(R^3))*(3*r(2)*r(3)/R^2);
% N(3,1) = (Ec.mu/(R^3))*(3*r(1)*r(3)/R^2);
% N(3,2) = (Ec.mu/(R^3))*(3*r(2)*r(3)/R^2);
% N(3,3) = (Ec.mu/(R^3))*(3*r(3)/R^2-1);

N(1,1) = (Const.MUe/(R^3))*(3*r(1)/R^2-1)+Ec.omega_e^2;
N(1,2) = (Const.MUe/(R^3))*(3*r(1)*r(2)/R^2);
N(1,3) = (Const.MUe/(R^3))*(3*r(1)*r(3)/R^2);
N(2,1) = (Const.MUe/(R^3))*(3*r(1)*r(2)/R^2);
N(2,2) = (Const.MUe/(R^3))*(3*r(2)/R^2-1)+Ec.omega_e^2;
N(2,3) = (Const.MUe/(R^3))*(3*r(2)*r(3)/R^2);
N(3,1) = (Const.MUe/(R^3))*(3*r(1)*r(3)/R^2);
N(3,2) = (Const.MUe/(R^3))*(3*r(2)*r(3)/R^2);
N(3,3) = (Const.MUe/(R^3))*(3*r(3)/R^2-1);

% Const.GM * Const.EM / R^3 * ((3/R^2 * r * r.') - eye(3));