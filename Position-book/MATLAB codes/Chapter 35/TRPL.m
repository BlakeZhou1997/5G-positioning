%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              TRPL.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lfs_dB,Lpe_dB,Lhh0,Lhh,LhhNG,Lvv0,Lvv,LvvNG]=TRPL(ff,epsilon,hh,lc,htt,hrr,dd)
% TRPL computes the near ground path loss given an infinitesimal dipole
% transmitter is radiating over a 1-D Gaussian dielectric rough surface.
%
%   [Lfs,Lpe,Lhh0,Lhh,LhhNG,Lvv0,Lvv,LvvNG]=TRPL(ff,epsilon,hh,lc,htt,hrr,dd)
%
%   INPUT:
%
%   ff=simulation frequency
%   epsilon=relative permittivity of lower medium
%   hh=rough surface rms height in lambda
%   lc=rough surface correlation length in lambda
%   htt=transmitter height in lambda
%   hrr=receiver height in lambda
%   dd=tranceiver distance in lambda
%
%   OUTPUT:
%
%   Lfs=free space path loss
%   Lpe=simplified plane earth path loss
%   Lhh0=path loss of the horizontally polarized field for horizontally
%   polarized incidence on a flat surface
%   Lhh=path loss of the horizontally polarized field for horizontally
%   polarized incidence
%   LhhNG=Lhh+knife-edge diffraction loss
%   Lvv0=path loss of the vertically polarized field for vertically
%   polarized incidence on a flat surface
%   Lvv=path loss of the vertically polarized field for vertically
%   polarized incidence
%   LvvNG=Lvv+knife-edge diffraction loss


c = 3e8;
lambda = c/ff;

ht=htt*lambda;
hr=hrr*lambda;
d=dd*lambda;

xt=0;
yt=0;
zt=ht;

xr=d;
yr=0;
zr=hr;

k=2*pi/lambda;
k1=k*sqrt(epsilon);

l=lc*lambda;
h=hh*lambda;

xd=xr-xt;
yd=yr-yt;

rhod=sqrt((xd)^2+(yd)^2);
phi_d=asin((yd)/rhod);

phi_sp=phi_d;
theta_sp=pi/2-atan((zt+zr)/rhod);

d_dir=sqrt(d^2+(ht-hr)^2);
d_ref=sqrt(d^2+(ht+hr)^2);
d_diff=d_ref-d_dir;
phase_diff=2*pi/lambda*d_diff;
Lfs=(2*k*d_dir)^2;
Lpe=(d^2/ht/hr)^2;

% Calculate the two ray path loss

Rhh0=Rh0(theta_sp,k,k1);
Rhh=Rh_eff(theta_sp,phi_sp,k,k1,h,l);
Rvv0=Rv0(theta_sp,k,k1);
Rvv=Rv_eff(theta_sp,phi_sp,k,k1,h,l);

Lhh0_excess=abs(1+d_dir/d_ref*Rhh0*exp(-1i*phase_diff))^-2;
Lhh_excess=abs(1+d_dir/d_ref*Rhh*exp(-1i*phase_diff))^-2;
Lvv0_excess=abs(1+d_dir/d_ref*Rvv0*exp(-1i*phase_diff))^-2;
Lvv_excess=abs(1+d_dir/d_ref*Rvv*exp(-1i*phase_diff))^-2;
L_ke=abs(0.5+0.8768*(ht+hr)/sqrt(lambda*d_dir))^-2;

Lfs_dB=10*log10(Lfs);
Lpe_dB=10*log10(Lpe);
Lhh0=10*log10(Lfs*max(Lhh0_excess,1));
Lhh=10*log10(Lfs*max(Lhh_excess,1));
LhhNG=10*log10(Lfs*max(Lhh_excess,1)*L_ke);
Lvv0=10*log10(Lfs*max(Lvv0_excess,1));
Lvv=10*log10(Lfs*max(Lvv_excess,1));
LvvNG=10*log10(Lfs*max(Lvv_excess,1)*L_ke);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         W_Spectrum_1D.m                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ww]=W_spectrum_1D(kx,l)
% define the spectrum of interest (Gaussian spectral density)

ww=l*sqrt(pi)*exp(-l^2*kx.^2/4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Rh0.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rh0]=Rh0(thetai,k,k1)
% 0th order solution: reflection coefficient for h-polarization

kz_sp=k*cos(thetai);
k1z_sp=sqrt(k1^2-(k*sin(thetai))^2);
Rh0=(kz_sp-k1z_sp)./(kz_sp+k1z_sp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Rv0.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rv0]=Rv0(thetai,k,k1)
% 0th order solution: reflection coefficient for v-polarization

kz_sp=k*cos(thetai);
k1z_sp=sqrt(k1^2-(k*sin(thetai))^2);
Rv0=(k1^2*kz_sp-k^2*k1z_sp)./(k1^2*kz_sp+k^2*k1z_sp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             omega_x.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=omega_x(thetai,phi_i,k,k1,h,l)
% calculate omega_x for a given spectrum

epsilon=(k1/k)^2;
kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

% integrate omega_x using numerical integration
n_kx=400;
aa=5/l;

kx=linspace(krhoi-aa,krhoi+aa,n_kx);
ky=kyi;
krho=sqrt(kx.^2+ky^2);
kz=sqrt(k^2-krho.^2);
k1z=sqrt(k1^2-krho.^2);
ww=W_spectrum_1D(krhoi-kx,l);
y11=k1z./(epsilon*kz+k1z).*(k1zi*kz+krhoi*kx);
y_ww=y11.*ww;
y_sum=trapz(kx,y_ww);
output=y_sum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             omega_y.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=omega_y(thetai,phi_i,k,k1,h,l)
% calculate omega_y for a given spectrum

epsilon=(k1/k)^2;
kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

% integrate omega_y using numerical integration
n_kx=400;
aa=5/l;

kx=linspace(krhoi-aa,krhoi+aa,n_kx);
ky=kyi;
krho=sqrt(kx.^2+ky^2);
kz=sqrt(k^2-krho.^2);
k1z=sqrt(k1^2-krho.^2);
ww=W_spectrum_1D(krhoi-kx,l);
y11=1./(kz+k1z);
y_ww=y11.*ww;
y_sum=trapz(kx,y_ww);
output=y_sum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             omega_z.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=omega_z(thetai,phi_i,k,k1,h,l)
% calculate omega_z for a given spectrum

epsilon=(k1/k)^2;
kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

% integrate omega_z using numerical integration
n_kx=400;
aa=5/l;

kx=linspace(krhoi-aa,krhoi+aa,n_kx);
ky=kyi;
krho=sqrt(kx.^2+ky^2);
kz=sqrt(k^2-krho.^2);
k1z=sqrt(k1^2-krho.^2);
ww=W_spectrum_1D(krhoi-kx,l);
y11=krho./(epsilon*kz+k1z).*(k1zi*kz+krhoi*kx);
y_ww=y11.*ww;
y_sum=trapz(kx,y_ww);
output=y_sum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Rh_eff.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Rh_eff(thetai,phi_i,k,k1,h,l)
% calculate Rh_eff, the effective reflection coefficient for a
% horizontally-polarized incident wave

epsilon=(k1/k)^2;
kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

y11=Rh0(thetai,k,k1);
y12=omega_y(thetai,phi_i,k,k1,h,l);
y13=h^2/2*k^2*(epsilon-1)*(1-y11^2);
y14=h^2/2*(k^2*(epsilon-1))^2*(1+y11)*y12/pi/(kzi+k1zi);
y_sum=y11+y13-y14;
output=y_sum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Rv_eff.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Rv_eff(thetai,phi_i,k,k1,h,l)
% calculate Rv_eff, the effective reflection coefficient for a
% vertically-polarized incident wave

epsilon=(k1/k)^2;
kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

y11=Rv0(thetai,k,k1);
y12=omega_x(thetai,phi_i,k,k1,h,l);
y13=omega_z(thetai,phi_i,k,k1,h,l);
y14=h^2/2*k^2*(epsilon-1)*(1-y11^2);
y15=h^2*(epsilon-1)*krhoi^2/(epsilon*kzi+k1zi)*((y11-1)*kzi-(y11+1)*k1zi);
y16=h^2/2*(epsilon-1)^2/(epsilon*kzi+k1zi)*((y11-1)*kzi*y12/pi+(y11+1)*krhoi*y13/pi);
y_sum=y11-y14-y15-y16;
output=y_sum;


