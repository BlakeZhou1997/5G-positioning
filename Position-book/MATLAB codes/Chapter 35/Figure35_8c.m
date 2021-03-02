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
% S=0.0502;
% C=0.4738;
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
rL=304.8;   %length of simulation area (along square side) in meters
Nn=100;   %number of nodes in the network
mc=1000;   %number of monte carlo iterations
n_mv=length(mv);

rng('shuffle');   %seeds the random number generator based on the current time so that RANDN produces different sequences of numbers

% Generate Nn*Nn uniformly distributed coordinates, find the pairwise
% distances and store them in d
for imc=1:mc
    d(imc,:,:)=dist(rL.*rand(2,Nn));
end

d=nonzeros(d);
d=reshape(d,[mc Nn,Nn-1]); %radial distance between transmitter and receiver in meters

% Calculate the geometric parameters in lambda to serve as inputs in TRPL
htt=ht/lambda;
hrr=hr/lambda;
dd=d/lambda;

% Calculate the break distance and the critical distance
d_b=sqrt((4*ht*hr/lambda-lambda/4)^2-(ht-hr)^2);
d_c=sqrt((12.5*ht*hr/lambda-lambda/12.5)^2-(ht-hr)^2);

for imv=1:n_mv
    for imc=1:mc
        for id1=1:Nn
            for id2=1:Nn-1

                % Calculate the path loss
                [Lfs(imv,imc,id1,id2),Lpe(imv,imc,id1,id2),Lhh0(imv,imc,id1,id2),Lhh(imv,imc,id1,id2),...
                LhhNG(imv,imc,id1,id2),Lvv0(imv,imc,id1,id2),Lvv(imv,imc,id1,id2),LvvNG(imv,imc,id1,id2)]=...
                TRPL(ff,epsilon(imv),hh,lc,htt,hrr,dd(imc,id1,id2));

                % Calculate the path loss for TE and TM polarization
                if (d(imc,id1,id2) <= d_b)
                    path_loss_Lhh0(imv,imc,id1,id2)=Lfs(imv,imc,id1,id2);
                    path_loss_LhhNG(imv,imc,id1,id2)=Lfs(imv,imc,id1,id2);
                    path_loss_Lvv0(imv,imc,id1,id2)=Lfs(imv,imc,id1,id2);
                    path_loss_LvvNG(imv,imc,id1,id2)=Lfs(imv,imc,id1,id2);
                elseif (d(imc,id1,id2) <= d_c)
                    path_loss_Lhh0(imv,imc,id1,id2)=Lhh0(imv,imc,id1,id2);
                    path_loss_LhhNG(imv,imc,id1,id2)=Lhh(imv,imc,id1,id2);
                    path_loss_Lvv0(imv,imc,id1,id2)=Lvv0(imv,imc,id1,id2);
                    path_loss_LvvNG(imv,imc,id1,id2)=Lvv(imv,imc,id1,id2);
                else
                    path_loss_Lhh0(imv,imc,id1,id2)=Lhh0(imv,imc,id1,id2);
                    path_loss_LhhNG(imv,imc,id1,id2)=LhhNG(imv,imc,id1,id2);
                    path_loss_Lvv0(imv,imc,id1,id2)=Lvv0(imv,imc,id1,id2);
                    path_loss_LvvNG(imv,imc,id1,id2)=LvvNG(imv,imc,id1,id2);
                end
            end
        end
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for imv=1:n_mv
    
    neighbors_Lhh0(imv)=length(find(path_loss_Lhh0(imv,:,:,:)<111))/mc/Nn;
    
    neighbors_Lvv0(imv)=length(find(path_loss_Lvv0(imv,:,:,:)<111))/mc/Nn;
    
    neighbors_LhhNG(imv)=length(find(path_loss_LhhNG(imv,:,:,:)<111))/mc/Nn;

    neighbors_LvvNG(imv)=length(find(path_loss_LvvNG(imv,:,:,:)<111))/mc/Nn;
    
end


figure (1)
plot(mv,neighbors_Lvv0);
hold on
plot(mv,neighbors_Lhh0);
plot(mv,neighbors_LvvNG);
plot(mv,neighbors_LhhNG);
legend('Flat ground, TM','Flat ground, TE','Rough ground, TM','Rough ground, TE');
xlabel('Volumetric moisture (cm^3/cm^3)');
ylabel('Average number of neighbors')
hold off


