%Shadowing variance
sigma1=[4,6,8];
%The variance of K is the sum of shadowing variance and penetration loss variance 
sigma0=sqrt(12.4^2+sigma1.^2);
%the range of K pdfs
x=-3*sigma0:0.1:3*sigma0;
%Initialize K pdf array for NLOS and LOS
kdb_nlos=zeros(length(sigma1),length(x));
kdb_los=kdb_nlos;

%Mean of K of LOS
u1=12;
%Mean of K of NLOS 
u0=u1-16.2;

%Compute pdfs of K
for index=1:length(sigma1)
kdb_nlos(index,:)=1/sqrt(2*pi)/sigma0(index)*exp(-(x-u0).^2/sigma0(index)^2);
kdb_los(index,:)=1/sqrt(2*pi)/sigma1(index)*exp(-(x-u1).^2/sigma1(index)^2);
end
figure
plot(x,kdb_nlos(1,:),'b--',x,kdb_nlos(2,:),'g-',x,kdb_nlos(3,:),'r-.',x,kdb_los(1,:),'b--',x,kdb_los(2,:),'g-',x,kdb_los(3,:),'r-','linewidth',2)
grid on
xlabel('K(dB)','fontsize',12)
legend(['\sigma_{sh}=',num2str(sigma1(1)),'dB'],['\sigma_{sh}=',num2str(sigma1(2)),'dB'],['\sigma_{sh}=',num2str(sigma1(3)),'dB']);
%Compute the threshold
kth=((sigma0.^2*u1-sigma1.^2*u0)-...
sqrt((sigma0.^2*u1-sigma1.^2*u0).^2-(sigma0.^2-sigma1.^2).*(sigma0.^2*u1^2-sigma1.^2*u0^2-sigma1.^2.*sigma0.^2.*log(sigma0./sigma1))))./(sigma0.^2-sigma1.^2)
%Compute the probability of detection 
PD=1-(0.5-0.5*erf((kth-u0)./sigma0/sqrt(2)))
%Compute the probability of false alarm
PF=1-(0.5-0.5*erf((kth-u1)./sigma1/sqrt(2)))

sh}=',num2str(sigma1(1)),'dB'],['\sigma_{sh}=',num2str(sigma1(2)),'dB'],['\sigma_{sh}=',num2str(sigma1(3)),'dB']);
%Compute the threshold
kth=((sigma0.^2*u1-sigma1.^2*u0)-...
sqrt((sigma0.^2*u1-sigma1.^2*u0).^2-(sigma0.^2-sigma