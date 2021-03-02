close all
clear all
clc

% Input arguments:

c = 3e8;  %speed of light in vacuum
ff=915e6;   %simulation frequency
lambda = c/ff;   %simuation wavelength in air
sigma=0.005;   %conductivity of the reflecting surface in mhos per meter
epsilon=15-60i*sigma*lambda;   %relative permittivity of lower medium (average ground)
hh=0.0113/lambda;   %rough surface rms height in lambda
lc=0.0739/lambda;   %rough surface correlation length in lambda
ht=0.05:0.05:0.95;   %transmitter height in meters
hr=ht;   %receiver height in meters
rL=304.8;   %length of simulation area (along square side) in meters
Nn=100;   %number of nodes in the network
mc=1000;   %number of monte carlo iterations
n_ht=length(ht);

rng('shuffle');   % seeds the random number generator based on the current time 
                  % so that RANDN produces different sequences of numbers

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

for iht=1:n_ht
    for imc=1:mc
        for id1=1:Nn
            for id2=1:Nn-1
        
                % Calculate the break distance and the critical distance
                d_b(iht)=sqrt((4*ht(iht)*hr(iht)/lambda-lambda/4)^2-(ht(iht)-hr(iht))^2);
                d_c(iht)=sqrt((12.5*ht(iht)*hr(iht)/lambda-lambda/12.5)^2-(ht(iht)-hr(iht))^2);

                % Calculate the path loss
                [Lfs(iht,imc,id1,id2),Lpe(iht,imc,id1,id2),Lhh0(iht,imc,id1,id2),Lhh(iht,imc,id1,id2),...
                LhhNG(iht,imc,id1,id2),Lvv0(iht,imc,id1,id2),Lvv(iht,imc,id1,id2),LvvNG(iht,imc,id1,id2)]=...
                TRPL(ff,epsilon,hh,lc,htt(iht),hrr(iht),dd(imc,id1,id2));

                % Calculate the path loss for TE and TM polarization
                path_loss_Lfs(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                path_loss_Lpe(iht,imc,id1,id2)=Lpe(iht,imc,id1,id2);
        
                if (d(imc,id1,id2) <= d_b(iht))
                    path_loss_Lhh0(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                    path_loss_Lhh(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                    path_loss_LhhNG(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                    path_loss_Lvv0(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                    path_loss_Lvv(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                    path_loss_LvvNG(iht,imc,id1,id2)=Lfs(iht,imc,id1,id2);
                elseif (d(imc,id1,id2) <= d_c(iht))
                    path_loss_Lhh0(iht,imc,id1,id2)=Lhh0(iht,imc,id1,id2);
                    path_loss_Lhh(iht,imc,id1,id2)=Lhh(iht,imc,id1,id2);
                    path_loss_LhhNG(iht,imc,id1,id2)=Lhh(iht,imc,id1,id2);
                    path_loss_Lvv0(iht,imc,id1,id2)=Lvv0(iht,imc,id1,id2);
                    path_loss_Lvv(iht,imc,id1,id2)=Lvv(iht,imc,id1,id2);
                    path_loss_LvvNG(iht,imc,id1,id2)=Lvv(iht,imc,id1,id2);
                else
                    path_loss_Lhh0(iht,imc,id1,id2)=Lhh0(iht,imc,id1,id2);
                    path_loss_Lhh(iht,imc,id1,id2)=Lhh(iht,imc,id1,id2);
                    path_loss_LhhNG(iht,imc,id1,id2)=LhhNG(iht,imc,id1,id2);
                    path_loss_Lvv0(iht,imc,id1,id2)=Lvv0(iht,imc,id1,id2);
                    path_loss_Lvv(iht,imc,id1,id2)=Lvv(iht,imc,id1,id2);
                    path_loss_LvvNG(iht,imc,id1,id2)=LvvNG(iht,imc,id1,id2);
                end
            end
        end
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for iht=1:n_ht
    
    neighbors_Lfs(iht)=length(find(path_loss_Lfs(iht,:,:,:)<111))/mc/Nn;
    
    neighbors_Lpe(iht)=length(find(path_loss_Lpe(iht,:,:,:)<111))/mc/Nn;
    
    neighbors_Lhh0(iht)=length(find(path_loss_Lhh0(iht,:,:,:)<111))/mc/Nn;
    
    neighbors_Lvv0(iht)=length(find(path_loss_Lvv0(iht,:,:,:)<111))/mc/Nn;

    neighbors_Lhh(iht)=length(find(path_loss_Lhh(iht,:,:,:)<111))/mc/Nn;

    neighbors_Lvv(iht)=length(find(path_loss_Lvv(iht,:,:,:)<111))/mc/Nn;
    
    neighbors_LhhNG(iht)=length(find(path_loss_LhhNG(iht,:,:,:)<111))/mc/Nn;

    neighbors_LvvNG(iht)=length(find(path_loss_LvvNG(iht,:,:,:)<111))/mc/Nn;
    
end


figure (1)
plot(ht,neighbors_Lfs);
hold on
plot(ht,neighbors_Lvv0);
plot(ht,neighbors_Lhh0);
plot(ht,neighbors_LvvNG);
plot(ht,neighbors_LhhNG);
legend('Free-space model','Two-ray model, TM','Two-ray model, TE','Proposed model, TM','Proposed model, TE');
xlabel('h_a (m)');
ylabel('Average number of neighbors')
hold off


