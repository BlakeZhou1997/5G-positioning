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
ht=0.05:0.1:1.95;   %transmitter height in meters
hr=ht;   %receiver height in meters
d=1:2000;   %radial distance between transmitter and receiver in meters
n_ff=length(ff);
n_d=length(d);
n_ht=length(ht);


for iff=1:n_ff
    for iht=1:n_ht
        for id=1:n_d
            
            % Calculate the geometric parameters in lambda to serve as inputs in TRPL
            htt(iff,iht)=ht(iht)/lambda(iff);
            hrr(iff,iht)=hr(iht)/lambda(iff);
            dd(iff,id)=d(id)/lambda(iff);
            
            % Calculate the break distance and the critical distance
            d_b(iff,iht)=sqrt((4*ht(iht)*hr(iht)/lambda(iff)-lambda(iff)/4)^2-(ht(iht)-hr(iht))^2);
            d_c(iff,iht)=sqrt((12.5*ht(iht)*hr(iht)/lambda(iff)-lambda(iff)/12.5)^2-(ht(iht)-hr(iht))^2);

            % Calculate the path loss
            [Lfs(iff,iht,id),Lpe(iff,iht,id),Lhh0(iff,iht,id),Lhh(iff,iht,id),...
            LhhNG(iff,iht,id),Lvv0(iff,iht,id),Lvv(iff,iht,id),LvvNG(iff,iht,id)]=...
            TRPL(ff(iff),epsilon(iff),hh(iff),lc(iff),htt(iff,iht),hrr(iff,iht),dd(iff,id));
            
            % Calculate the path loss for TE and TM polarization
            if (d(id) <= d_b(iff,iht))
                path_loss_hh(iff,iht,id)=Lfs(iff,iht,id);
                path_loss_vv(iff,iht,id)=Lfs(iff,iht,id);
            elseif (d(id) <= d_c(iff,iht))
                path_loss_hh(iff,iht,id)=Lhh(iff,iht,id);
                path_loss_vv(iff,iht,id)=Lvv(iff,iht,id);
            else
                path_loss_hh(iff,iht,id)=LhhNG(iff,iht,id);
                path_loss_vv(iff,iht,id)=LvvNG(iff,iht,id);
            end
        end   
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for iff=1:n_ff
    for iht=1:n_ht

    iLhh(iff,iht)=find(squeeze(path_loss_hh(iff,iht,:))>111,1);
    coverage_range_hh(iff,iht)=d(iLhh(iff,iht));
    
    iLvv(iff,iht)=find(path_loss_vv(iff,iht,:)>111,1);
    coverage_range_vv(iff,iht)=d(iLvv(iff,iht));
    
    end
end


figure (1)
semilogy(ht,coverage_range_vv(1,:));
hold on
semilogy(ht,coverage_range_hh(1,:));
semilogy(ht,coverage_range_vv(2,:));
semilogy(ht,coverage_range_hh(2,:));
semilogy(ht,coverage_range_vv(3,:));
semilogy(ht,coverage_range_hh(3,:));
legend('f = 315 MHZ, TM','f = 315 MHZ, TE','f = 915 MHZ, TM','f = 915 MHZ, TE','f = 2.4 GHZ, TM','f = 2.4 GHZ, TE');
xlabel('h_a (m)');
ylabel('Maximum link range (m)')
hold off



