% -------------------------------------------------------------------------
% Chapter 23:
% An overview on Global Positioning Techniques for Harsh Environments
% Example 2: Correlation between two spreading codes with noise,
%            exploiting integration time extension
% Author: Nicola Linty, Politecnico di Torino
% Based on a code by Marco Pini, Istituto Superiore Mario Boella (Chapter 21)
% -------------------------------------------------------------------------

%--- clean up the environment
clear
close all
clc

%--- Variables definition
Rc = 0.5e6;                              % Code Rate (chip/s)
Fs = 8e6;                                % Sampling Frequency (Hz)

%--- define PRN codes
cLc = [1 -1 1 -1 -1 -1 -1 1 -1 ...
         1 -1 -1 -1 1 -1 -1 1 1 1 -1];

L = length(cLc);                         % Code Length (chip)
N = floor(Fs*L/Rc);                      % Code Length (samples)

%--- Sample the local code
k = 0:N-1;                               
cLocSampled = cLc(floor(k*Rc/Fs)+1); % cLoc code sampled at Fs
%--- Number of sums for high sensitivity
Nsums = 10;
cLocSampled = repmat(cLocSampled, [1 Nsums]);

% - Generate the incoming code with M=20 periods of cLoc, sample and shift
M = 20;
cIn = repmat(cLc, [1 M]);
k = 0:M*N-1;
cInSampled = cIn(floor(k*Rc/Fs)+1);   % cIn code sampled @ Fs

%--- Samples of the incoming code with a code-phase shift
Delay = 4*Fs/Rc;                         % Code Delay (samples)
cInSampledShift = circshift(cInSampled, [0 Delay]);

%--- Add AWGN noise
sigmaAWGN = 1;
% sigmaAWGN = 12.5;
cInSampledNoise = cInSampledShift + sigmaAWGN * randn(1, M*N);
% -------------------------------------------------------------------------

%--- Correlate the two sequences of samples
CorrFull = zeros(1, Nsums*N); % initialize the variable
for index = 0:Nsums*N-1
   %--- correlate the codes
   CorrFull(index+1) = cInSampledNoise(1+index:2*N+index) * cLocSampled(1:2*N)';                    
end
%--- sum the correlation results
Corr = sum(reshape(CorrFull, N, Nsums)');

%--- Plot correlation functions
xAxis = [0:(N-1)] ./ Fs .* Rc;            % Prepare x-axis (chip)

figure
plot(xAxis, Corr, '.-k')
grid on
xlabel('Delay (chip)')
ylabel('Correlation')
title('Code cross-correlation')
axis tight
