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
ht=0.05:0.1:1.95;   %transmitter height in meters
hr=ht;   %receiver height in meters
d=1:2000;   %radial distance between transmitter and receiver in meters
n_d=length(d);
n_ht=length(ht);

% Calculate the geometric parameters in lambda to serve as inputs in TRPL
htt=ht/lambda;
hrr=hr/lambda;
dd=d/lambda;

for iht=1:n_ht
    for id=1:n_d
        
        % Calculate the break distance and the critical distance
        d_b(iht)=sqrt((4*ht(iht)*hr(iht)/lambda-lambda/4)^2-(ht(iht)-hr(iht))^2);
        d_c(iht)=sqrt((12.5*ht(iht)*hr(iht)/lambda-lambda/12.5)^2-(ht(iht)-hr(iht))^2);

        % Calculate the path loss
        [Lfs(iht,id),Lpe(iht,id),Lhh0(iht,id),Lhh(iht,id),...
        LhhNG(iht,id),Lvv0(iht,id),Lvv(iht,id),LvvNG(iht,id)]=...
        TRPL(ff,epsilon,hh,lc,htt(iht),hrr(iht),dd(id));

        % Calculate the path loss for TE and TM polarization
        path_loss_Lfs(iht,id)=Lfs(iht,id);
        path_loss_Lpe(iht,id)=Lpe(iht,id);
        
        if (d(id) <= d_b(iht))
            path_loss_Lhh0(iht,id)=Lfs(iht,id);
            path_loss_Lhh(iht,id)=Lfs(iht,id);
            path_loss_LhhNG(iht,id)=Lfs(iht,id);
            path_loss_Lvv0(iht,id)=Lfs(iht,id);
            path_loss_Lvv(iht,id)=Lfs(iht,id);
            path_loss_LvvNG(iht,id)=Lfs(iht,id);
        elseif (d(id) <= d_c(iht))
            path_loss_Lhh0(iht,id)=Lhh0(iht,id);
            path_loss_Lhh(iht,id)=Lhh(iht,id);
            path_loss_LhhNG(iht,id)=Lhh(iht,id);
            path_loss_Lvv0(iht,id)=Lvv0(iht,id);
            path_loss_Lvv(iht,id)=Lvv(iht,id);
            path_loss_LvvNG(iht,id)=Lvv(iht,id);
        else
            path_loss_Lhh0(iht,id)=Lhh0(iht,id);
            path_loss_Lhh(iht,id)=Lhh(iht,id);
            path_loss_LhhNG(iht,id)=LhhNG(iht,id);
            path_loss_Lvv0(iht,id)=Lvv0(iht,id);
            path_loss_Lvv(iht,id)=Lvv(iht,id);
            path_loss_LvvNG(iht,id)=LvvNG(iht,id);
            end
    end
end

% Find the maximum coverage range for a transmitted power of 10 dBm and a sensitivity of -101 dBm using various models    
for iht=1:n_ht
    
    iLfs(iht)=find(path_loss_Lfs(iht,:)>111,1);
    coverage_range_Lfs(iht)=d(iLfs(iht));
    
    iLpe(iht)=find(path_loss_Lpe(iht,:)>111,1);
    coverage_range_Lpe(iht)=d(iLpe(iht));

    iLhh0(iht)=find(path_loss_Lhh0(iht,:)>111,1);
    coverage_range_Lhh0(iht)=d(iLhh0(iht));
    
    iLvv0(iht)=find(path_loss_Lvv0(iht,:)>111,1);
    coverage_range_Lvv0(iht)=d(iLvv0(iht));

    iLhh(iht)=find(path_loss_Lhh(iht,:)>111,1);
    coverage_range_Lhh(iht)=d(iLhh(iht));
    
    iLvv(iht)=find(path_loss_Lvv(iht,:)>111,1);
    coverage_range_Lvv(iht)=d(iLvv(iht));

    iLhhNG(iht)=find(path_loss_LhhNG(iht,:)>111,1);
    coverage_range_LhhNG(iht)=d(iLhhNG(iht));
    
    iLvvNG(iht)=find(path_loss_LvvNG(iht,:)>111,1);
    coverage_range_LvvNG(iht)=d(iLvvNG(iht));
    
end


figure (1)
semilogy(ht,coverage_range_Lfs);
semilogy(ht,coverage_range_Lvv0);
hold on
semilogy(ht,coverage_range_Lhh0);
semilogy(ht,coverage_range_LvvNG);
semilogy(ht,coverage_range_LhhNG);
legend('Two-ray model, TM','Two-ray model, TE','Proposed model, TM','Proposed model, TE');
xlabel('h_a (m)');
ylabel('Maximum link range (m)')
hold off



