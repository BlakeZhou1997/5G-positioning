% -------------------------------------------------------------------------
% Chapter 21 Example 1: Correlation between two spreading codes with noise
% -------------------------------------------------------------------------

Rc = 0.5e6;                              % Code Rate [chip/s]
Fs = 8e6;                                % Sampling Frequency [Hz]

c_loc = [1 -1 1 -1 -1 -1 -1 1 -1 ...
         1 -1 -1 -1 1 -1 -1 1 1 1 -1];   % PRN code

% sample the spreading codes
L = length(c_loc);                       % Code Length [chip]
N = floor(Fs*L/Rc);                      % Code Length [samples]

% ------------------------- Sample the local code -------------------------
k = 0:N-1;
c_loc_sampled = c_loc(floor(k*Rc/Fs)+1); % c_loc code sampled @ Fs
% -------------------------------------------------------------------------

% - Generate the incoming code with 3 periods of c_loc, sample and shift --
c_in = [c_loc c_loc c_loc]; k = 0:3*N-1;
c_in_sampled = c_in(floor(k*Rc/Fs)+1);   % c_in code sampled @ Fs


D = 4*Fs/Rc;                            % Code Delay [samples]

% samples of the incoming code with a code-phase shift
c_in_sampled_shift = circshift(c_in_sampled,[0 D]);
% -------------------------------------------------------------------------

% correlate the two sequences of samples
for index=0:N-1,
    % correlate the codes...
    Corr(index+1) = c_in_sampled_shift(1+index:N+index)*c_loc_sampled(1:N)';
end

% plot correlation functions
x_axis = [0:(N-1)]./Fs.*Rc;            % Prepare x-axis [chip]

figure(1),plot(x_axis,Corr,'.-k'), grid on;
xlabel('Delay [chip]')
ylabel('Correlation')
title('PRN code correlation')
