function Ec = Earth_constants
% USAGE:
%   Ec = Earth_constants
%
% OUTPUT:
%   Ec = structure with the following fields:
%       omega_e: Earth rate magnitude (Sideral Rate) (rad/s) (Ref[2])
%       M: Earth mass (kg) (Ref[2])
%       G: Universal Gravitational constant (Ref[2])
%       mu: G*M
%       R0: Earth Equatorial Radius (m) (Ref[2])
%       e: Earth ellipticity (aka flattening) (Ref[2])
%       J2,J3: Earth figure constants (Ref[1]).
%
% NOTES:
% Ellipticity is equal to flattening, f, and related to eccentricity via:
%   eccentricity^2 = f(2-f).
%
% REFERENCES:
%   [1] P.G.Savage, "Strapdown Analytics, Vols 1&2", 2000.
%   [2] www.wikipedia.org
%   [3] Britting,K.R.,"Inertial Navigation System Analysis", John Wiley &
%   Sons, New York, 1971
%   [4] Brookes,C.J., "Evaluation of Odd Zonal Harmonics in the Earth's
%   Gravitational Potential", Proc.Math.Phys.Sci.,V446,No1926,1994,149-168

% Earth rate magnitude (Sideral Rate) (Ref[2]):
Ec.omega_e = 1e0*(2*pi)/((23*60+56)*60+4.091);

% Earth mass (Ref[2]):
Ec.M = 5.9742e24; % Kg

% Universal Gravitational constant (Ref[2]):
Ec.G = 6.6742e-11; % N.m^2.kg^-2

Ec.mu = Ec.G*Ec.M;

% Earth Equatorial Radius (Ref[2]):
Ec.R0 = 6.378137e6; % m %/
Ec.SM_AXIS = 6.378137e6; % m %//!Semi-major axis of Ellipsoid in m

Ec.SMIN_AXIS = 6356752.3142; % m %//!Semi-major axis of Ellipsoid in m

% Earth ellipticity:
% NB: Equal to flattening, f, and related to eccentricity via:
%   eccentricity^2 = f(2-f).
Ec.e = 1/298.257223563; % (Ref[1]) (it is called also FLATTENING)

Ec.WIE_E=7292115e-11;        %earth's rotation rate

Ec.NORMAL_GRV=9.7803253359;  %earth's NORMAL GRAVITY  (different from  Standard GRAVITY) //!Theoretical (Normal) Gravity at the equator (on Ellipsoid)

Ec.STD_GRV=9.80665;  %earth's STANDARD GRAVITY 

Ec.GRV_CONS=0.00193185265241; %Normal Gravity Constant 

Ec.M_FAKTOR=0.00344978650684; %%//!combined term for gravity calculation

Ec.E_SQR=Ec.e*(2-Ec.e);%=0.00669437999014; %%!First Eccentricity squared

Ec.L1              = 1575420000;             % GPS L1 Frequency  

Ec.C               = 299792458.0;            % Speed of light                                        	

Ec.lambda=0.1902936728; % for L1

Ec.F               = -4.442807633e-10;       % Relativistic correction term constant     

Ec.WE              = 7.292115e-5;            % WGS84 Earth rotation rate    

Ec.WEDOT           =  7.292115147e-5;        % WGS84 Value of earth's rotation rate                  


% Earth gravitational harmonic coefficients (Ref[1]).
Ec.J2 = 1.082627e-3; % Ref[1]
Ec.J3 = -2.5327e-6; % Ref[1]
Ec.J4 = -1.8e-6; % Ref[3]
Ec.J5 = -0.247e-6;% Ref[4]