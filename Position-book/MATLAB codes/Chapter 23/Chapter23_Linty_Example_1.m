% -------------------------------------------------------------------------
% Chapter 23:
% An overview on Global Positioning Techniques for Harsh Environments
% Example 1: Correlation between two spreading codes with noise,
%            effect of a longer code
% Author: Nicola Linty, Politecnico di Torino
% Based on a code by Marco Pini, Istituto Superiore Mario Boella (Chapter 21)
% -------------------------------------------------------------------------

%--- Clean up the environment
clear
close all
clc

%--- Variables definition
Rc = 0.5e6;                              % Code Rate (chip/s)
Fs = 8e6;                                % Sampling Frequency (Hz)

%--- define PRN codes
cLoc1 = [1 -1 1 -1 -1 -1 -1 1 -1 ...
         1 -1 -1 -1 1 -1 -1 1 1 1 -1];   
cLoc2 = [ 1 1 -1 1 -1 1 1 -1 -1 1 ...
         1 1 1 1 -1 -1 -1 1 1 -1 1 ...
         1 1 -1 1 -1 1 -1 -1 -1 -1 ];

% cLoc = cLoc1;
cLoc = cLoc2;

L = length(cLoc);                        % Code Length (chip)
N = floor(Fs*L/Rc);                      % Code Length (samples)

%--- Sample the local code
k = 0:N-1;                               
cLocSampled = cLoc(floor(k*Rc/Fs)+1); % cLoc code sampled at Fs

%--- Generate the incoming code with 3 periods of cLoc, sample and shift
cIn = [cLoc cLoc cLoc];
k = 0:3*N-1;
cInSampled = cIn( floor( k*Rc/Fs ) + 1 );  % cIn code sampled at Fs

%--- Samples of the incoming code with a code-phase shift
Delay = 4*Fs/Rc;                        % code Delay (samples)
cInSampledShift = circshift(cInSampled, [0 Delay]);

%--- Add AWGN noise
sigmaAWGN = 1;
% sigmaAWGN = 4;
% sigmaAWGN = 8;

cInSampledNoise = cInSampledShift + sigmaAWGN * randn(1, 3*N);

%--- Correlate the two sequences of samples
Corr = zeros(1, N); % initialize the variable
for index = 0:N-1
   %--- correlate the codes
   Corr(index+1) = cInSampledNoise(1+index:N+index) * cLocSampled(1:N)';                    
end    

%--- Plot correlation functions
xAxis = [0:(N-1)] ./ Fs .* Rc;         % Prepare x-axis (chip)

figure,
plot(xAxis, Corr, '.-k'),
grid on
xlabel('Delay (chip)')
ylabel('Correlation')
title('PRN code correlation')
axis tight