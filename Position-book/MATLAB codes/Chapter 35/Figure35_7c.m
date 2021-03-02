close all
clear all
clc

% Input arguments:

c = 3e8;  %speed of light in vacuum
ff=915e6;   %simulation frequency
lambda = c/ff;   %simuation wavelength in air
mv=0:0.02:0.50;   %volumetric moisture content
%selected soil type is silt loam (sand(%)=17.16,silt(%)=63.84,clay(%)=19)
S=0.1716;
C=0.19;
epsilon0=8.854e-12;
rho_s=2.66;   %N.R. Peplinski, F.T. Ulaby, M.C. Dobson, "Dielectric properties of soils in the 0.3-1.3-GHz range"
rho_b=1.37;   %http://www.pedosphere.com/resources/bulkdensity/worktable_us.cfm
epsilon_s=(1.01+0.44*rho_s)^2-0.062;
beta_p=1.2748-0.519*S-0.152*C;
beta_dp=1.33797-0.603*S-0.166*C;
alpha=0.65;
sigma_eff=0.0467+0.2204*rho_b-0.4111*S+0.6614*C;
epsilon_fw_p=4.9+(80.1-4.9)./(1+(0.58e-10*ff).^2);
epsilon_fw_dp=(0.58e-10*ff*(80.1-4.9))./(1+(0.58e-10*ff).^2)+(sigma_eff*(rho_s-rho_b))./(2*pi*ff*epsilon0*rho_s*(mv+1e-15));
epsilon=1.15.*(1+rho_b/rho_s.*(epsilon_s^alpha-1)+mv.^beta_p.*epsilon_fw_p.^alpha-mv).^(1/alpha)-0.68-1i.*(mv.^beta_dp.*epsilon_fw_dp.^alpha).^(1/alpha);   %relative permittivity of lower medium
hh=0.0113/lambda;   %rough surface rms height in lambda
lc=0.0739/lambda;   %rough surface correlation length in lambda
ht=0.10;   %transmitter height in meters
hr=ht;   %receiver height in meters
d=1:500;   %radial distance between transmitter and receiver in meters
n_mv=length(mv);
n_d=length(d);

% Calculate the geometric parameters in lambda to serve as inputs in TRPL
htt=ht/lambda;
hrr=hr/lambda;
dd=d/lambda;

% Calculate the break distance and the critical distance
d_b=sqrt((4*ht*hr/lambda-lambda/4)^2-(ht-hr)^2);
d_c=sqrt((12.5*ht*hr/lambda-lambda/12.5)^2-(ht-hr)^2);

for imv=1:n_mv
    for id=1:n_d
        
        % Calculate the path loss
        [Lfs(imv,id),Lpe(imv,id),Lhh0(imv,id),Lhh(imv,id),...
        LhhNG(imv,id),Lvv0(imv,id),Lvv(imv,id),LvvNG(imv,id)]=...
        TRPL(ff,epsilon(imv),hh,lc,htt,hrr,dd(id));

        % Calculate the path loss for TE and TM polarization
            if (d(id) <= d_b)
                path_loss_hh0(imv,id)=Lfs(imv,id);
                path_loss_hhNG(imv,id)=Lfs(imv,id);
                path_loss_vv0(imv,id)=Lfs(imv,id);
                path_loss_vvNG(imv,id)=Lfs(imv,id);
            elseif (d(id) <= d_c)
                path_loss_hh0(imv,id)=Lhh0(imv,id);
                path_loss_hhNG(imv,id)=Lhh(imv,id);
                path_loss_vv0(imv,id)=Lvv0(imv,id);
                path_loss_vvNG(imv,id)=Lvv(imv,id);
            else
                path_loss_hh0(imv,id)=Lhh0(imv,id);
                path_loss_hhNG(imv,id)=LhhNG(imv,id);
                path_loss_vv0(imv,id)=Lvv0(imv,id);
                path_loss_vvNG(imv,id)=LvvNG(imv,id);
            end
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for imv=1:n_mv
    
    iLhh0(imv)=find(path_loss_hh0(imv,:)>111,1);
    coverage_range_hh0(imv)=d(iLhh0(imv));
    
    iLvv0(imv)=find(path_loss_vv0(imv,:)>111,1);
    coverage_range_vv0(imv)=d(iLvv0(imv));

    iLhhNG(imv)=find(path_loss_hhNG(imv,:)>111,1);
    coverage_range_hhNG(imv)=d(iLhhNG(imv));
    
    iLvvNG(imv)=find(path_loss_vvNG(imv,:)>111,1);
    coverage_range_vvNG(imv)=d(iLvvNG(imv));
    
end


figure (1)
semilogy(mv,coverage_range_vv0);
hold on
semilogy(mv,coverage_range_hh0);
semilogy(mv,coverage_range_vvNG);
semilogy(mv,coverage_range_hhNG);
legend('Flat ground, TM','Flat ground, TE','Rough ground, TM','Rough ground, TE');
xlabel('Volumetric moisture (cm^3/cm^3)');
ylabel('Maximum link range (m)')
hold off


