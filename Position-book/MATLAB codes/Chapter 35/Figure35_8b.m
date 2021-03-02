close all
clear all
clc

% Input arguments:

c = 3e8;  %speed of light in vacuum
ff=[315e6 915e6 2400e6];   %simulation frequency
lambda = c./ff;   %simuation wavelength in air
sigma=0.005;   %conductivity of the reflecting surface in mhos per meter
epsilon=15-60i*sigma*lambda;   %relative permittivity of lower medium (average ground)
hh=0.0113./lambda;   %rough surface rms height in lambda
lc=0.0739./lambda;   %rough surface correlation length in lambda
ht=0.05:0.05:0.95;   %transmitter height in meters
hr=ht;   %receiver height in meters
rL=304.8;   %length of simulation area (along square side) in meters
Nn=100;   %number of nodes in the network
mc=1000;   %number of monte carlo iterations
n_ff=length(ff);
n_ht=length(ht);

rng('shuffle');   % seeds the random number generator based on the current 
                  % time so that RANDN produces different sequences of numbers

% Generate Nn*Nn uniformly distributed coordinates, find the pairwise
% distances and store them in d
for imc=1:mc
    d(imc,:,:)=dist(rL.*rand(2,Nn));
end

d=nonzeros(d);
d=reshape(d,[mc Nn,Nn-1]); %radial distance between transmitter and receiver in meters

for iff=1:n_ff
    for iht=1:n_ht
        for imc=1:mc
            for id1=1:Nn
                for id2=1:Nn-1
        
                    % Calculate the geometric parameters in lambda to serve as inputs in TRPL
                    htt(iff,iht)=ht(iht)/lambda(iff);
                    hrr(iff,iht)=hr(iht)/lambda(iff);
                    dd(iff,imc,id1,id2)=d(imc,id1,id2)/lambda(iff);
                
                    % Calculate the break distance and the critical distance
                    d_b(iff,iht)=sqrt((4*ht(iht)*hr(iht)/lambda(iff)-lambda(iff)/4)^2-(ht(iht)-hr(iht))^2);
                    d_c(iff,iht)=sqrt((12.5*ht(iht)*hr(iht)/lambda(iff)-lambda(iff)/12.5)^2-(ht(iht)-hr(iht))^2);

                    % Calculate the path loss
                    [Lfs(iff,iht,imc,id1,id2),Lpe(iff,iht,imc,id1,id2),Lhh0(iff,iht,imc,id1,id2),Lhh(iff,iht,imc,id1,id2),...
                    LhhNG(iff,iht,imc,id1,id2),Lvv0(iff,iht,imc,id1,id2),Lvv(iff,iht,imc,id1,id2),LvvNG(iff,iht,imc,id1,id2)]=...
                    TRPL(ff(iff),epsilon(iff),hh(iff),lc(iff),htt(iff,iht),hrr(iff,iht),dd(iff,imc,id1,id2));

                    % Calculate the path loss for TE and TM polarization
                    if (d(imc,id1,id2) <= d_b(iff,iht))
                        path_loss_hh(iff,iht,imc,id1,id2)=Lfs(iff,iht,imc,id1,id2);
                        path_loss_vv(iff,iht,imc,id1,id2)=Lfs(iff,iht,imc,id1,id2);
                    elseif (d(imc,id1,id2) <= d_c(iff,iht))
                        path_loss_hh(iff,iht,imc,id1,id2)=Lhh(iff,iht,imc,id1,id2);
                        path_loss_vv(iff,iht,imc,id1,id2)=Lvv(iff,iht,imc,id1,id2);
                    else
                        path_loss_hh(iff,iht,imc,id1,id2)=LhhNG(iff,iht,imc,id1,id2);
                        path_loss_vv(iff,iht,imc,id1,id2)=LvvNG(iff,iht,imc,id1,id2);
                    end
                end
            end
        end
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for iff=1:n_ff
    for iht=1:n_ht
    
        neighbors_hh(iff,iht)=length(find(path_loss_hh(iff,iht,:,:,:)<111))/mc/Nn;
        neighbors_vv(iff,iht)=length(find(path_loss_vv(iff,iht,:,:,:)<111))/mc/Nn;
    
    end
end

figure (1)
plot(ht,neighbors_vv(1,:));
hold on
plot(ht,neighbors_hh(1,:));
plot(ht,neighbors_vv(2,:));
plot(ht,neighbors_hh(2,:));
plot(ht,neighbors_vv(3,:));
plot(ht,neighbors_hh(3,:));
legend('f = 315 MHZ, TM','f = 315 MHZ, TE','f = 915 MHZ, TM','f = 915 MHZ, TE','f = 2.4 GHZ, TM','f = 2.4 GHZ, TE');
xlabel('h_a (m)');
ylabel('Average number of neighbors')
hold off


