% Matlab codes to generate PDP for S-V channel

Lmbd = 1/200; % Cluster arrival rate
lmbd = 1/5; % Ray arrival rate

Gmm = 60; % Cluster decay rate
gmm = 20; % Ray decy rate

Th = 300; % Maximum channel delay

T(1) = 0; % Cluster arrival time 'T'
num_cluster = 1; % Number of clusters

while 1        
    % Generate next cluster arrival time
    temp_T = T(num_cluster) + exprnd(1/Lmbd); 
    
    % Maximum channel delay < Th
    if temp_T < Th
        num_cluster = num_cluster + 1;
        T(num_cluster) = temp_T;
    else
        break;
    end
end

PDP(1) = 1; % Power delay profile 'PDP'
num_MPC = 0; % Number of multipath components 

for m = 1:num_cluster    
    num_MPC = num_MPC + 1;    
    t(num_MPC) = T(m); % Multipath delay 't'
    % Generate PDP 
    PDP(num_MPC) = PDP(1)*exp(-T(m)/Gmm)*...
        exp(-(t(num_MPC)-T(m))/gmm);
    
    while 1          
        % Ray arrival time inside clusters
        temp_t = t(num_MPC) + exprnd(1/lmbd);
        
        % Maximum channel delay less than Th
        if temp_t < Th
            num_MPC = num_MPC + 1;            
            t(num_MPC) = temp_t; % Save multipath delay  
            % Save PDP
            PDP(num_MPC) = PDP(1)*exp(-T(m)/Gmm)*...
                exp(-(t(num_MPC)-T(m))/gmm);
        else
            break;            
        end
    end    
end

% Normalized PDP
PDP = PDP/sum(PDP);

stem(t, PDP);
xlabel('Delay (ns)');
ylabel('Normalized PDP for S-V channel');

