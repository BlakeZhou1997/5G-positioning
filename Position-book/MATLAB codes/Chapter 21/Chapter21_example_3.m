% -------------------------------------------------------------------------
% Chapter 21 Example 3: Serial search acquisition algorithm.
% -------------------------------------------------------------------------

clear all; close all; clc

Rc = 0.5e6;                              % Code Rate [chip/s]
Fs = 8e6;                                % Sampling Frequency [Hz]
Fc = 2e6;                                % Carrier Frequency [Hz]

c_loc = [1 -1 1 -1 -1 -1 -1 1 -1 ...
         1 -1 -1 -1 1 -1 -1 1 1 1 -1];   % PRN code

L = length(c_loc);                       % Code Length [chip]
N = floor(Fs*L/Rc);                      % Code Length [samples]


% --------- Generate the incoming signal, with 3 periods of c_loc ---------
c_in = [c_loc c_loc c_loc];

% sample the spreading codes
k = 0:3*N-1;
c_in_sampled = c_in(floor(k*Rc/Fs)+1);   % c_in code sampled @ Fs

D = 4*Fs/Rc;                             % code delay [samples]
sigma = 0.5;

% modulate the spreading code sampled @ Fs and add noise
c_in_sampled_noise = circshift(c_in_sampled,[0 D]).*cos(2*pi*Fc.*k/Fs) ...
    + (sigma)*randn(1,3*N);
% -------------------------------------------------------------------------

% ----------------------- Generate the local code -------------------------
k = 0:N-1;                               % samples vector
c_loc_sampled = c_loc(floor(k*Rc/Fs)+1); % c_loc code sampled @ Fs
% -------------------------------------------------------------------------

for idx_freq = 1:10,

    fd = Fc-50e3 + 10e3*idx_freq;    % local carrier frequency [Hz]
    loc_cos = cos(2*pi*fd.*k/Fs);    % In-phase local carrier
    loc_sin = sin(2*pi*fd.*k/Fs);    % Quadrature local carrier

    % prepare the local signal, multipling local code and carriers
    Loc_cos = c_loc_sampled.*loc_cos;
    Loc_sin = c_loc_sampled.*loc_sin;

    for idx_shift = 0:N-1,
        % correlate incoming and local signals
        Corr(idx_freq,idx_shift+1) = ...
            (c_in_sampled_noise(1+idx_shift:N+idx_shift)*Loc_cos').^2 + ...
            (c_in_sampled_noise(1+idx_shift:N+idx_shift)*Loc_sin').^2;
    end
end

% plot correlation functions
x_axis = [0:(N-1)]./Fs.*Rc;              % Prepare x-axis [chip]
y_axis = Fc-5e3 + 1e3.*[1:10]

figure(1), surf(x_axis,y_axis,abs(Corr)), shading interp,
xlabel('Delay [chip]')
ylabel('Freq [MHz]')
title('Cross Ambiguity Function')

