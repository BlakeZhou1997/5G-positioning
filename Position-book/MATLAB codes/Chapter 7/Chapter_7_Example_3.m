close all, clear all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EXAMPLE 7.3 :Simple correlation based TOA estimation method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses Signal Processing Toolbox
%------------------------------------------------------------------------
% the first arriving peak is estimated by comparing the output of 
% cross-correlator module with a certain threshold.
% the threshold can be set based on the measured power of noise.
                     % calculate the cross correlation between received signal and transmit signal
%------------------------------------------------------------------------
% to evaluate this method, we use examples 1 & 2 to generate data
%------------------------------------------------------------------------
% generate multipath wireless channel
% base band model of Saleh-Valenzuela (S-V) channel model
% The model is implemented as a tap delay line structure 
% The model assumes that multipath components arrive in clusters
% the cluster arrival rate is described by a Poisson process 
% The average power of subsequent clusters is assumed to decay exponentially. 
% Each cluster is composed of many multipath rays, 
% whose arrival times are also described by Poisson process
% the average ray power in any given cluster is assumed to decay exponentially 
% the phase of each ray is uniformly distributed over [0, 2*pi].
ts = 50 ;                                        % sampling time (nsec) (BW = 1/ts)
t0=550.4;                                     % First cluster arrival occurs 
fc0=1e-3 * [5800];                       % first center frequency per GHz
fc= [fc0 fc0+1e9/ts fc0+2e9/ts];  % center frequency per GHz, fc must be in GHz in accordance with ts
channel_length=100;                   % number of samples in channel model
SNR_dB=10;                               % Signal To Noise Ratio in dB
%------------------------------------------------------------------------
%S-V model parameters
Gam = 550;                             % Cluster decay factor
gamma = 670;                         % Ray decay factor
st1 = .4;                                   % Standard deviation of log-normal variable for cluster fading
st2 = .4;                                   % Standard deviation of log-normal  variable for ray fading
mu_const = 0.3;
cluster_nom=10;                     % number of cluster in generated channel  
ray_per_cluster_nom=5;         % number of rays in each cluster
max_channel_length=4000;    % maximum channel duration in nSec
max_cluster_length=400;        % maximum cluster duration in nSec
%------------------------------------------------------------------------
% get one realizations of continus time channel impulse response  
cluter_vect=(max_channel_length-t0)*rand(1,cluster_nom);
cluster_matrix=repmat(cluter_vect,ray_per_cluster_nom,1);
ray_matrix=max_cluster_length*rand(ray_per_cluster_nom,cluster_nom);
ln_xi_vect=st1*randn(size(cluter_vect));      % set cluster fading 
ln_xi_matrix=repmat(ln_xi_vect,ray_per_cluster_nom,1);
ln_beta_matrix = st2*randn(size(ray_matrix)); % set ray fading 
final_delays=ray_matrix+cluster_matrix;       % time of arrival matrix
mu_matrix = cluster_matrix/Gam-ray_matrix/gamma - mu_const; 
                          % ray decay because of path
                          % loss ( function of cluster delay and ray decay)
h_val_matrix = exp(-1*(ln_xi_matrix+ln_beta_matrix+mu_matrix));
                          % generate ray coeeficent
temp_t=final_delays(:);
temp_h=h_val_matrix(:);
np = cluster_nom*ray_per_cluster_nom;    % number of rays (or paths) for this realization
[sort_temp_t,sort_ix] = sort(temp_t);    % sort in ascending time order
t_ct = sort_temp_t;
h = temp_h(sort_ix);
                 % multiply each tap with a uniform random phase [0,2*pi]
h_ct = h .*  exp(-sqrt(-1)* 2 * pi * rand(size(h)));
% Outputs
%   h_ct is returned  a random realization of the channel impulse response 
%   t_ct is organized as h_ct, but holds the time instances (in nsec) of 
%   the paths whose signed amplitudes are stored in h_ct
%   t0 is the arrival time of the first cluster 
%   np is the number of paths 
%------------------------------------------------------------------------
% convert continuous-time channel model h_ct to N-times oversampled samples
% ts is the desired time resolution
% hN will be produced with time resolution ts / N.
% It is up to the user to then apply any filtering and/or complex 
% downconversion and then decimate by N to finally obtain 
% an channel impulse response (CIR) at time resolution ts.
% N*fs = N/ts is the intermediate sampling frequency before decimation
N =128;                         % make N a power of 2 to facilitate efficient multi-stage decimation
Nfs = N / ts;                   % bandwidths of channel (GHz)    
t_max = max(t_ct(:));           % maximum time value across all channels
h_len = 1 + floor(t_max * Nfs); % number of time samples at resolution ts/N
hN = zeros(h_len,1);
np_k = np;                      % number of paths in this channel
t_Nfs = 1 + floor(t_ct(1:np) * Nfs);  % vector of quantized time 
for n = 1:np
    hN(t_Nfs(n),1) = hN(t_Nfs(n),1) + h_ct(n,1);
end
%------------------------------------------------------------------------
% To obtain our desired complex base band equivalent channel model sampled,
% first shift the whole spectrum by the desired center frequency,
% apply low pass filtering (with cut off frequency BW/2 MHz) 
% and then decimate down by N to obtain our desired complex base band equivalent channel model
% repeat this process for each subband
for k=1:length(fc),
  hN_b_shift(:,k)=hN .* exp(-sqrt(-1)*2*pi* fc(k)*(0:length(hN)-1)*ts/N).';
        % Resample and correct for 1/N scaling imposed by decimation
  h_b(:,k) = N * resample(hN_b_shift(:,k), 1, N);
end
                % Restrict the length of channel to a given channel_length
h_b =[h_b ; zeros(channel_length,size(h_b,2))];
h_b=h_b(1:channel_length,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we use an OFDM based system with these parameters
%      bandwidth         : 1/ts (Based on Nyquist rate),
%      subcarrier number : sample_nom,
%      Modulation        : QPSK modulation with random data
% For some high resolution approaches, multiple symbols are generated.
% the distance between two consecutive symbols is selected properly to 
% avoid Inter Symbol Interference (ISI) due to the multipath effect. 
% Using IFFT module, we convert these data to time domain. 
sample_nom=64;                  % number of sample per symbol
symbol_nom=120;                 % number of symbol 
%------------------------------------------------------------------------
%  generate Transmit signal
%  OFDM symbol with QPSK modulation
                               % generate QPSK subcarrier in frequency domain
sf=ones(sample_nom,1).* ...                
           exp(sqrt(-1)*((pi*floor(4*rand(sample_nom,1))/2)+pi/4));
tx_symbol_freq_domain=[sf(sample_nom/2+1:sample_nom);
    sf(1:sample_nom/2)];
s=ifft(tx_symbol_freq_domain,sample_nom);  % convert signal to time domain
tx_matrix=repmat(tx_symbol_freq_domain,1,symbol_nom);
tx_sig=real([s;zeros(127,1);s;zeros(127,1)]);
%------------------------------------------------------------------------
% applying channel to Transmit signal
% The multipath channel is applied by convolving the transmitted signal 
% with the complex base band channel model generated in example 7.1.
% the additive noise is modeled with a white complex Gaussian random variable. 
% The power of this noise is extracted from the power 
% of received signal and the desired SNR.
for k=1:length(fc),   % applying channel by  convolving tx signal with CIR
    single_h_b=[h_b(:,k); zeros(128-size(h_b,1),1)];
    multiband_channel_out(:,k)=conv(s,single_h_b);
end
                      % measure the power of signal to caculate the power of noise from SNR
rx_sig_power=sum(abs(multiband_channel_out(:)).^2) ...
               /numel(multiband_channel_out(:));
                      % calculate the power of noise
component_noise_sigma=sqrt(2)/2*10^(-SNR_dB/20)*sqrt(rx_sig_power);
                      % generate number of symbol as a packet of data
                      % here we assume there is no Inter Symbol Interference (ISI) exists in the system 
channel_out_matrix=repmat(multiband_channel_out(:,1),1,symbol_nom);
                      % add white quassian additive noise to signal
rx_sig_matrix = channel_out_matrix + ...
    component_noise_sigma*(randn(size(channel_out_matrix))+ ...
    sqrt(-1)*randn(size(channel_out_matrix)));
          % in multiband scenario we assume we send eachsymbol in one subband defined in channel model
multiband_rx_sig = multiband_channel_out + ...
    component_noise_sigma*(randn(size(multiband_channel_out))+ ...
    sqrt(-1)*randn(size(multiband_channel_out)));
 
rx_sig=real([rx_sig_matrix(:,1);rx_sig_matrix(:,2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcor_output=xcorr(rx_sig_matrix(:,1),s);
xcor_output=xcor_output(length(rx_sig_matrix(:,1))+1:end);
                 % Compare the absolute value of xcor output with a cetain threshold to detect first received ray
                     % this threshold can be selected based on the noise  power
ray_index_xcor=find((abs(xcor_output) ) > ...
                         4*component_noise_sigma);
if ~isempty(ray_index_xcor)
   xcor_estimated_TOA=ray_index_xcor(1)*ts
else
   xcor_estimated_TOA=NaN
end
figure, plot(ts*(1:length(xcor_output)),abs(xcor_output),'k' );
hold on;
plot(ts*(1:length(xcor_output)), ...
               component_noise_sigma*10*ones(size(xcor_output)),'k-.');
xlabel('nSec');
title('cross correltor output');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
