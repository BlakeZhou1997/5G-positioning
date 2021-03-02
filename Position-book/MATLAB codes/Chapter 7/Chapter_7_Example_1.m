clear all, close all, clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 7.1  Channel Model Generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses Signal Processing Toolbox
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
figure, stem(t_ct,abs(h_ct),'k');
% xlabel('nSec'), xlim([0,4000]);
title('continus time channel impulse response');
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
figure, plot(ts*(1:channel_length),abs(h_b(:,1)),'k');
xlabel('nSec');
title('channel impulse response for BW=20MHz');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
