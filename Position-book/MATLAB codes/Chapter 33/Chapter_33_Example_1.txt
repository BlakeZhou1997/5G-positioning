%Robert MacCurdy - rbm7@cornell.edu - September 2010
%Example 1: Numerical simulation of a time of arrival receiver

%This example uses Matlab to generate a carrier signal, modulates that
%carrier signal with a Gold code, adds Gaussian noise, down converts the
%signal into I/Q baseband signals, cross-correlates the signals with the
%Gold code template, and uses a threshold detector to indicate signal
%detection. The code illustrates the importance of I/Q baseband processing
%in a communications system that is not phase synchronous, and provides a
%way of investigating the impact of clock frequency offsets in the tag and
%receiver.

clear %MUST DO THIS for proper operation!

%%%%%%%%%%%%%USER-DEFINED CONSTANTS%%%%%%%%%%%

CarrierFreq=150.0e6;    %(Hz): tag carrier frequency (enter desired ideal carr freq - offset is handled below)
ChipRate=1e6;           %(Hz): pn sequence chip rate
SNR=20;                 %(dB) Signal to Noise Ratio
MixerFreq=150.0e6;      %(Hz): LO mixer freq
DeltaFreq=0;            %(Hz): the amount by which the tag's carrier freq deviates
                        %from the LO freq.
PhaseDiff=(0*pi)/8;     %(Radians): phase shift of the carrier relative to the LO
LpfStopFreq=2e6;        %(Hz): Cutoff frequency of the image rejection low pass 
                        %filter after the mixer. (note the way this is implemented prevents aliasing, even if it is desired)
BBsr=2.8169e6;          %(Hz): Base-band sample rate. Simulates the A/D sampling
                        %occurs occurs after the mixer (equivalent to 
                        %resampling in systems that digitize at IF and 
                        %perform digital downconversion)
                        
SimSampRt=1e9;          %(Hz): The simulation sample rate. This number must be at
                        %least twice as high as the highest frequency in
                        %the simulation in order to abide by the Nyquist
                        %criterion. Higher sample rates make the plots look
                        %better                       

%randn('state',0)       %this line can be enabled to always use the same
                        %random numbers in generating the contaminating 
                        %noise that is added to the signal

ThreshScale=1.4;        %(Dimensionless) Used to set the threshold for our detector (see 
                        %below).

%%%%%%NO USER-MODIFIABLE CONSTANTS BELOW%%%%%%%%

sim('Chapter_33_Example_1_pn_codes'); %Create the PN Gold code "Gold0" used below. 
%NOTE: if the "Chapter_33_Example_1_pn_codes.mdl" Simulink file is not available, it can be
%recreated in the following way: Run Simulink by typing "simulink" in the 
%command window. Create a new empty, simulink model by selecting new->model
%from the Simulink file menu. Insert a "Gold Sequence Generator" block, 
%which is available in the "Communications Blockset->Comm Sources->Sequence
%Generators" library. Also add a "To Workspace" block, available in the
%"Simulink->Sinks" library. Connect these two blocks with a signal. Double
%click on the To Workspace block and change the Variable Name to "Gold0".
%Change the "Save Format:" parameter to "Array".
%Click Apply and then Ok. Double click the Gold Sequence Generator block
%and change the Preferred polynomial (1) to [11 2 0], the Initial states
%(1) to 1, the Preferred Polynomial (2) to [11 8 5 2 0], the Initial States
%(2) to 1, the Sequence Index to 8, and the Shift to 0. Click the 
%Simulation->Configuration Parameters... menu and set the Start time to 0,
%the Stop time to 1456, the Solver "Type" to "Fixed-Step", and the "Solver"
%to "discrete (no continuous states)". Click Apply and Ok. 

numcyc = CarrierFreq/ChipRate;  %number of cycles per chip

%generate carrier (possibly w/ freq offset) and Lo
Tau1 = 1/(CarrierFreq);       %wavelength of desired carrier
t=0:(1/SimSampRt):length(Gold0)*numcyc*Tau1-(1/SimSampRt);
carr=sin((2*pi*CarrierFreq-DeltaFreq)*t+PhaseDiff);
LoI=sin(2*pi*MixerFreq*t);
LoQ=cos(2*pi*MixerFreq*t);

%scale and upsample the pn sequence appropriately so we can easily multiply
%the signals
pn=round(Gold0-.5);    %want values of +/- 1
pn = resample(pn',length(t),length(pn));

%mix the pn sequence with the carrier, as the tag would do with BPSK
carr=carr.*pn;

%downconvert quadrature signals
I=carr.*LoI;
Q=carr.*LoQ;

%downsample to Baseband sample rate (note that decimate includes a filter stage
%to avoid aliasing)
BBI=decimate(I,round(SimSampRt/BBsr));
BBQ=decimate(Q,round(SimSampRt/BBsr));
BBpn=decimate(pn,round(SimSampRt/BBsr));

%filter signal - this simulates the action of any potential post-mixer
%image blocking LP filter (note that the previous decmiate operation already
%did this for us, but don't have that in an analog system)
LPfilt_ratio=LpfStopFreq/(BBsr/2);  %do some range checking to prevent the filter function from failing
if (LPfilt_ratio>=1); LPfilt_ratio=.99; end
if (LPfilt_ratio<=0); LPfilt_ratio=.01; end
[B,A] = butter(8,LPfilt_ratio,'low');   %get coefficients for a Butterworth filter
BBI=filtfilt(B,A,BBI);  %filtfilt runs the filter forward and backward to eliminate filter phase delay
BBQ=filtfilt(B,A,BBQ);

%introduce noise to the baseband signal (to simulate corruption with noise 
%at all stages in signal chain)
rmssig=norm(BBI+BBQ.*1i)/sqrt(length(BBI));  %measure the RMS value of the signal
rmsnoise=rmssig/(10^(SNR/20));              %compute desired RMS value of noise
noise=rmsnoise*randn(1,1*length(BBI));      %produce the noise

%Add the noise to the baseband signal.
BBIn=noise+BBI;
BBQn=noise+BBQ;
tBB=linspace(0,(length(BBIn)-1)/BBsr,length(BBIn));

%At this point, we have our downconverted, sampled Baseband signal that is
%contaminated with noise. We can look at it:
figure(1);
plot(tBB,BBIn,'b',tBB,BBQn,'g');

%Look at the cross-correlation of pn w I&Q individually. Notice that the
%cross-correlation peak for each channel depends on the PhaseDiff
%parameter even when the SNR is high? Not good. Imagine if we only looked
%at the I channel, which is all we'd have wthout the quadrature signal in
%the mixer. It would dropout each time the carrier and LO were pi/2 out of
%phase!
[Xcorr_PN_I,Lags_PN_I] = xcorr([BBpn,BBpn,BBpn],BBIn);
[Xcorr_PN_Q,Lags_PN_Q] = xcorr([BBpn,BBpn,BBpn],BBQn);
figure(2);
%We hav to select a smaller range of the full cross-correlation to plot,
%since xcorr performs a full cross-correlation, not a circular
%cross-correlation (it pads with zeros). This causes the ends of the
%cross-correlation output to look lower in magnitude than they would be if
%they were part of a continuous signal, which is what we want to show.
strt=find(Lags_PN_I==0)+1;
endpt=strt+2*(length(BBpn)-1);
subplot1=subplot(2,1,1);plot(Lags_PN_I(strt:endpt),20*log10(Xcorr_PN_I(strt:endpt)),'b');
subplot2=subplot(2,1,2);plot(Lags_PN_Q(strt:endpt),20*log10(Xcorr_PN_Q(strt:endpt)),'g');
ylim(subplot1,[-20 80]);
ylim(subplot2,[-20 80]);

%now look at the cross-correlation of the full complex baseband signal
%with the PRN template. Notice that this correlation peak is independent of
%phase. 
[XcorrBB,LagsBB] = xcorr(([BBpn+BBpn*1i,BBpn+BBpn*1i,BBpn+BBpn*1i]),(BBIn+BBQn*1i));
strt=find(LagsBB==0)+1;
endpt=strt+2*(length(BBpn)-1);
figure3=figure(3);
clf(figure3);   %clear the figure, since we're using "hold all" to handle 
                %multiple plots
% Create axes
axes1 = axes('Parent',figure3);
hold(axes1,'all');
plot(axes1,LagsBB(strt:endpt),20*log10(abs(XcorrBB(strt:endpt))));

%finally, let's implement a simple threshold detector. we'll use the median
%value of the cross-correlation output, and require that the correlation
%peak exceed it by a certain margin for a detection to be registered.
Threshold=ThreshScale*20*log10(median(abs(XcorrBB(strt:endpt))));
% Create line indicating the threshold
plot(axes1,[LagsBB(strt),LagsBB(endpt)],[Threshold,Threshold],'r','LineWidth', 2);

