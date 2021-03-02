%This algorithm uses the phase difference of two signals received by an
%antenna array to estimate the rician K factor.
clear;
%parameters of the K estimator
a=0.0738; b=0.4842; c=0.1256; beta=1.004; gamma=-0.6121;
%number of samples
N=200;
%Each K is estimated for 100 times to obtain the root mean squre error of
%the estimation.
iteration=100;

%The ture values of K
K=0:50;
%initialize the array for estimated K
Kphi=zeros(iteration,length(K));
%r1 and r2 are signals received from antenna 1 and antenna 2
%Phase of r1's LOS component
phi_los1=2*pi*rand;
%LOS componenet of r1
rlos1=exp(1j*phi_los1);
%phase shift between the LOS componets of r1 and r2 
phase_shift=2*pi*rand;
%LOS component of r2
rlos2=exp(1j*(phi_los1+phase_shift));

for iteration_index=1:iteration
    %NLOS components of r1 and r2
    rnlos1=1/sqrt(2)*(randn(N,1)+randn(N,1)*1j);
    rnlos2=1/sqrt(2)*(randn(N,1)+randn(N,1)*1j);
    %initialize phase difference variance
    sigma2=zeros(size(K));
    
    for index=1:length(K)
        %Generate r1 and r2
        r1=sqrt(K(index))*rlos1+rnlos1;
        r2=sqrt(K(index))*rlos2+rnlos2;
        
        %phi1 and phi2 are wrapped in the range [-pi,pi]
        phi1=angle(r1);
        phi2=angle(r2);
                
        %phi1_2pi and phi2_2pi are wrapped in the range [0, 2pi]
        phi1_2pi=wrapto2pi(phi1);
        phi2_2pi=wrapto2pi(phi2);
        
        %Compare the variance of phi1, phi2 and phi1_2pi, phi2_2pi, 
        %and choose the wrapping with a smaller variance
        if var(phi1)>var(phi1_2pi)
            phi1_new=phi1_2pi;
        else
            phi1_new=phi1;
        end
        if var(phi2)>var(phi2_2pi)
            phi2_new=phi2_2pi;
        else
            phi2_new=phi2;
        end
        
        %Compute the phase difference variance
        sigma2(index)=var(phi2_new-phi1_new);
    end
    %Compute K from phase difference variance
    for index=1:length(K)
        if sigma2(index)>0.232
            Kphi(iteration_index,index)=(-b+sqrt(b^2-4*a*(c-1/sigma2(index))))/(2*a);
        else
            Kphi(iteration_index,index)=(1/sigma2(index)-gamma)/beta;          
        end
    end 
end

%Compute the root mean square error of K estimation
rms_phi=sqrt(mean((Kphi-ones(iteration,1)*K).^2));
%Compute the bias of K estimation
bias_phi=mean((Kphi-ones(iteration,1)*K));

figure
subplot(2,1,1)
plot(K,rms_phi,'k*')
xlabel('K')
ylabel('Root mean square error')
title('Root mean squre error of K estimation based on phase difference')
grid on
subplot(2,1,2)
plot(K,bias_phi,'o')
xlabel('K')
ylabel('Bias')
title('Bias of K estimation based phase difference')
grid on
d
    end 
end

%Compute the root mean square error of K estimation
rms_phi=sqrt(mean((Kphi-ones(iteration,1)*K).^2));
%Compute the bias of K estimation
bias_phi=mean((Kphi-ones(iteration,1)*K));

figure
subplot(2,1,1)
plot(K,rms_phi,'k*')
xlabel('K')
ylabel('Root mean square error')
title('Root mean squre error of K estimation