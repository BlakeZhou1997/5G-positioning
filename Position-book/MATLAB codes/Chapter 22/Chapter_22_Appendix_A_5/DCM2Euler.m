function Euler = DCM2Euler(R)
% DESCRIPTION:
%   Returns Euler angles from b to n frame DCM
%
% USAGE:
%   Euler = DCM2Euler(R)
%
% INPUT:
%   R:  b to n frame DCM
%
% OUTPUT:
%   Euler:  structure with the following fields:
%     azimuth:    vehicle azimuth, measured clockwise from North, in rad.
%     pitch:      vehicle elevation, measured upward from horizontal, in rad.
%     roll:       vehicle roll, in radians

Euler.yaw = atan2(R(2,1),R(1,1)); % eq. 3-7
Euler.pitch = -asin(R(3,1)); % eq. 3-6
Euler.roll = atan2(R(3,2),R(3,3)); % eq. 3-5

% End of function: DCM2Euler
return
