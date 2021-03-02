%---------------------------------------------------------------
% Chapter 21 Example 4: FFT-based CAF evaluation
%---------------------------------------------------------------

clear all
close all
clc

%------------------- Load of gold code (GPS)

load('GPS_code.mat')
disp('Contents of workspace after loading file:')
whos

%------------------- parameters of the carrier
fd=-2750;                        % Doppler frequency
FIF=4e6;                         % Intermediate frequency
BW=2*FIF;
phi0=0;
FRF=1575.42e6;                   % Radio frequency (GPS L1)
%------------------- parameters of the PRN code
Nchip=1023;                      % number of chips
Tcode=1e-3;                      % Code period
Tchip_nom=Tcode/Nchip;
Tchip=Tchip_nom/(1+fd/FRF);    % chip period with Doppler effect

%------------------- parameters related to noise and SIS power
CN0dBHz=48;                    % carrier to noise ratio
CN0=10^(CN0dBHz/10);
SNR=CN0/BW;                    % Signal to noise ratio
sigmanoise=1;
sigmanoise2=sigmanoise^2;
ampSIS=sqrt(2)*sqrt(SNR*sigmanoise2);

%------------------- sampling frequency
fs=2.1*BW;                     % sampling frequency
Ts=1/fs;

%------------------- simulation paramters
N_Tcode=3;
T_fin=Tchip*Nchip*N_Tcode;     % Duration of the simulated signal
time=0:Ts:T_fin;
Ntime=length(time);

%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% Signal generation

disp('Signal generation')

% N_Tcode code periods are generated

code_rep=code1023;
for indc=1:N_Tcode-1
    code_rep=[code_rep, code1023];
end

time_floor=floor(time/Tchip);

icounter=1;
for index=1:Ntime-1
    codesignal(index)=code_rep(icounter);
    if time_floor(index+1) > time_floor(index)
        icounter=icounter+1;
    end
end
% codesignal is the PRN code (with N_Tcode code periods)

carrier=cos(2*pi*(FIF+fd)*time(1:end-1)+phi0);
SIS_clean=ampSIS*codesignal.*carrier;
SIS=SIS_clean+sigmanoise*randn(1,length(carrier));
% SIS is the received signal with noise (with no data)

% ------------- Extraction of a SIS segment of 1ms
t1num=300;
t1=t1num*Tchip;         % Simulation of an arbitrary delay
t2=t1+Tcode;
delay=Tcode-t1;

N1=floor(t1/Ts);
N2=floor(t2/Ts);

SIS_seg=SIS(N1:N2);
NSIS=length(SIS_seg);
% SIS_seg is a segment of the received signal (1 code period)

%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% CAF evaluation (FFT in the time domain)

disp('CAF evaluation (FFT in the time domain)')
fmax=5000;                         % maximum Doppler frequency
fdbar=-fmax:300:+fmax;             % Frequency bins for CAF
LF=length(fdbar);
timeCAF_bin=0:Ts:Tcode;            % Time bins for CAF

DFTc=conj(fft(codesignal(1:length(timeCAF_bin))));

for indf=1:LF
    eloc=exp(1i*2*pi*(FIF+fdbar(indf))*timeCAF_bin);
    qSe=SIS_seg.*eloc;
    DFTq=fft(qSe);
    CAF_in_TD(:,indf)=ifft(DFTc.*DFTq);
end
% CAF_in_TD is the CAF evaluated by FFT in the time domain

% --- CAF sections
[m1,asc1]=max(abs(CAF_in_TD));
[m2,asc2]=max(m1);
CAFtime=CAF_in_TD(:,asc2);
[m3,asc3]=max(CAFtime);
CAFfreq=CAF_in_TD(asc3,:);

%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% CAF evaluation (FFT in the frequencye domain)

disp('CAF evaluation (FFT in the frequencye domain)')

% --------- Definition of the prefiter parameters
BW_FIR=4*fmax;                    % value arbitrarily chosen
T_FIR=1/BW_FIR;
N_dec=round(T_FIR/Ts);
bFIR=ones(1,N_dec)/N_dec;

% --------- Parameters for signal decimation
Num_col=floor(NSIS/N_dec);
Ntot=Num_col*N_dec;
ND=round(N_dec/2);

% --------- Mixer with complex exponential (for down-conversion)
SIS_with_mixer=SIS_seg.*exp(-1i*2*pi*FIF*time(1:length(SIS_seg)));
num_delay=length(SIS_seg);

perc0=10;
Dec_segments=floor(num_delay/perc0);
perc=perc0;

for indd=1:num_delay
    disp_num=mod(indd,Dec_segments);
    if disp_num==0
        disp([num2str(perc),'% of evaluated cells'])
        perc=perc+perc0;
    end
    code_del=codesignal(1+indd-1:length(timeCAF_bin)+indd-1);
    SIS_car=code_del.*SIS_with_mixer;
    SIS_car_FIR=filter(bFIR,1,SIS_car);
    % down-conversion
    SIS_car_FIR_matrix=reshape(SIS_car_FIR(1:Ntot),[N_dec,Num_col]);
    SIS_dec=SIS_car_FIR_matrix(ND,:);
    % decimation (digital sampling)
    LS=length(SIS_dec);
    CAF_in_FD(indd,:)=fft([SIS_dec zeros(1,LS)]);
    % zero padding arbitrarily chosen (can be modified)
end
% CAF_in_FD is the CAF evaluated by FFT in the frequency domain

%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% Figures

% CAF in the time domain
mesh(fdbar,timeCAF_bin,abs(CAF_in_TD))
title('CAF in the time domain')
xlabel('Doppler frequency (Hz)')
ylabel('Code delay (s)')

set(gcf, 'PaperPosition', [0 0 18 12]);
set(gcf, 'PaperSize', [18 12]);
saveas(gcf, 'CAF_TD', 'fig') %Save figure
saveas(gcf, 'CAF_TD', 'eps') %Save figure

figure
subplot(211)
plot(timeCAF_bin,abs(CAFtime))
legend(['True delay =',num2str(delay),' s'])
grid on
subplot(212)
plot(fdbar,abs(CAFfreq))
legend(['True Doppler freq. =',num2str(fd),' Hz'])
grid on

% CAF in the frequency domain
figure
[nrow,ncol]=size(CAF_in_FD);
ncdiv2=round(ncol/2);
cc1=CAF_in_FD(:,1:ncdiv2);
cc2=CAF_in_FD(:,ncdiv2+1:end);
CAF_in_FD_sim=[cc2 cc1];    % matrix reordering

Tdur=Tcode+LS*Ts*N_dec;
delf_fft=1/Tdur;
CAFfreq_FD=0:1:ncol-1;
CAFfreq_FD=CAFfreq_FD-ncdiv2;
CAFfreq_FD=CAFfreq_FD*delf_fft;

CAF_in_FD_sim_ud=flipud(CAF_in_FD_sim); % matrix reordering
mesh(CAFfreq_FD,timeCAF_bin,abs(CAF_in_FD_sim_ud))
title('CAF in the frequency domain')
xlabel('Doppler frequency (Hz)')
ylabel('Code delay (s)')

set(gcf, 'PaperPosition', [0 0 18 12]);
set(gcf, 'PaperSize', [18 12]);
saveas(gcf, 'CAF_FD', 'fig') %Save figure
saveas(gcf, 'CAF_FD', 'eps') %Save figure

[m1,asc1]=max(abs(CAF_in_FD_sim_ud));
[m2,asc2]=max(m1);
CAFtime=CAF_in_FD_sim_ud(:,asc2);
[m3,asc3]=max(CAFtime);
CAFfreq=CAF_in_FD_sim_ud(asc3,:);

figure
subplot(211)
plot(timeCAF_bin,abs(CAFtime))
legend(['True delay =',num2str(delay),' s'])
grid on
subplot(212)
plot(CAFfreq_FD,abs(CAFfreq))
legend(['True Doppler freq. =',num2str(fd),' Hz'])
grid on
