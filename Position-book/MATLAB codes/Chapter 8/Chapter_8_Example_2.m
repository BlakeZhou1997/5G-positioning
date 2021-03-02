% clear screen
clc;
% clear worksapce
clear;

% L: number of channel taps
L = 12;
% Lmt: number of purse noise samples
Lmt = 5;

% Td: channel delay spread
Td = 70;

% N: Simulation number
N = 20000;

% Nb: number of subbands
Nb = 6;

% m: Nakagami-m parameter
m = 1;

% snr: signal to noise ratio 
snr = 0:2:30;
SNR = 10.^((snr/10));

% t: tap delays of the discrete time channel 
t = linspace(0, Td, L);

% decay_factor: last tap energy / first tap energy
decay_factor = 100.0;

% factor_a: PDP fading factor
factor_a = log(decay_factor)/Td;

% hpw: exponential PDP
hpw = exp(-factor_a*t);

hwind = zeros(2*Lmt+1);

Pe = zeros(size(snr));
Pe_ub = zeros(L, length(snr));

% channel paths
h = zeros(Nb,L+2*Lmt);

% n0: noise samples
n0 = zeros(Nb,L+2*Lmt);

for ln = 1:N
    
    % print the progress of simulation
    if mod(ln, 1000) == 0
        ln
    end


    % lb: index of subband
    for lb = 1:Nb        
        % thta: channel tap phases
        thta = 2*pi*rand(1,L);

        % h: Nakagami-m channel taps
        A = m;
        B = hpw/m;
        dlt = gamrnd(A, B, 1, L);

        h(lb,:) = [zeros(1,Lmt), sqrt(dlt).*exp(sqrt(-1)*thta), zeros(1,Lmt)];
    end

    % noise with unity variance
    n0 = sqrt(1/2)*randn(size(h))+sqrt(-1/2)*randn(size(h));

    % lsnr: index of snr
    for lsnr = 1:length(snr)
        
        % add noise to channel taps
        h1 = h*sqrt(SNR(lsnr)) + n0;

        h1 = h1.*conj(h1);

        % noncoherent combining of subbands
        h2 = sum(h1, 1);

        % calculate energy for L samples
        for lw = 1:length(h1)-L+1
            hwind(lw) = sum(h2(lw:lw+L-1));
        end
        
        % union bound of the mistiming probability
        for lw = 1:Lmt
            if (hwind(lw)>hwind(Lmt+1))
                Pe_ub(lw, lsnr) = Pe_ub(lw, lsnr) + 1;
            end
            
            if (hwind(lw+Lmt+1)>hwind(Lmt+1))
                Pe_ub(lw, lsnr) = Pe_ub(lw, lsnr) + 1;
            end
        end       
        

        % I: estimated TOA
        [C,I] = max(hwind);

        % Pe: number of timing errors
        if (I~=Lmt+1)
            Pe(lsnr) = Pe(lsnr) + 1;
        end
    end

end

% output result
figure;
semilogy(snr, Pe/N);
hold on;



