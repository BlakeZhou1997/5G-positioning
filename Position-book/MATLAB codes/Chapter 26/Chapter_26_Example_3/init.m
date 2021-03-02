%%%%%%%% Initialization %%%%%%%%%%%%%%%

%%% Parameters (can be modified):
placement='dete';       % 'rand' (random), or 'dete' (deterministic) mode
two_hop='y';            % 'y' - include 2-hop neighbors, 'n' - don't include
area=0:0.2:20;          % range of the deployment area
sigma=0.02*(max(area)-min(area)); % std. deviation of the distance measurements

%%% Node placement (complex notation: x=(a,b)=a+j*b):
if placement=='rand' % expect 4 anchors
    Niter=4;  % number of iterations (minimum 2)
    Na=4;     % number of anchors (minimum 4)
    Nn=25+Na; % number of nodes (minimum 5)
    Np=100;   % number of particles (preferable 100-1000)
    R=0.40*(max(area)-min(area));       % transmission radius
    % First 4 anchors (close to 4 edges):
    x(1)=1+j*1;                         % anchor 
    x(2)=1.5+j*0.96*max(area);          % anchor
    x(3)=0.95*max(area)+j*0.5;          % anchor
    x(4)=0.93*max(area)+j*0.91*max(area); % anchor
    % Rest of the nodes (random):
    for n=5:Nn 
        x(n)=((max(area)-min(area))*rand(1,1)+min(area))+j*((max(area)-min(area))*rand(1,1)+min(area)); 
    end
elseif placement=='dete' % determenistic placement (including anchors)
    Niter=2; 	% number of iterations (minimum 2)
    Na=3;       % number of anchors
    Nn=2+Na;    % number of nodes
    Np=1000;    % number of particles
    R=0.4*(max(area)-min(area));       % transmission radius
    % Set x manually (first Na values are anchors):
	x(1)=4+j*13;        % anchor
	x(2)=4+j*6;         % anchor
  	x(3)=17+j*14;       % anchor
   	x(4)=6.5+j*8;       % target
   	x(5)=14.5+j*8;      % target
end

%%% Define distance matrix:
for m=1:Nn
    for n=1:Nn
        if m~=n & abs(x(m)-x(n))<=R 
            dist(m,n)=abs(x(m)-x(n));
        else
            dist(m,n)=-1; % no detection
        end
    end
end

%%% No-link status (only 2-hop neighbors) for bounding box:
for m=1:Nn
    for n=1:Nn
        if dist(m,n)<0 % no detection
            for k=Na+1:Nn
                if dist(m,k)>0 & dist(n,k)>0 no_link_status(m,n)=1; end
            end
        else no_link_status(m,n)=0; end
    end
end

%%% Define bounded box for particles (include 1-hop and 2-hop anchors):
for n=Na+1:Nn
	x_min(n)=-10*max(area)-j*10*max(area); % any large negative number
    x_max(n)=10*max(area)+j*10*max(area);  % any large positive number
    kk=0;
    for m=1:Na
        if dist(m,n)>0 % direct distance to anchor node
           dist2anchor(n,m)=dist(m,n); 
        elseif no_link_status(m,n)==1 % 2-hop distance to anchor node
           d_min=10*max(area); % any large positive number
           for k=(Na+1):Nn                  
                if d_min>(dist(m,k)+dist(k,n)) & dist(m,k)>0 & dist(k,n)>0
                    % we need upper bound, the sum
                    dist2anchor(n,m)=dist(m,k)+dist(k,n); d_min=dist2anchor(n,m); 
                end
           end
           if d_min==10*max(area) dist2anchor(n,m)=-2; end
        else dist2anchor(n,m)=-2; % no 1 and 2 hop connections
        end
        % Find bounds using min-max (use 3*sigma to make sure that the true position is within the bounds)
        if dist2anchor(n,m)>0
           x_min(n)=max(real(x(m))-dist2anchor(n,m)-3*sigma,real(x_min(n)))...
           +j*max(imag(x(m))-dist2anchor(n,m)-3*sigma,imag(x_min(n)));
           x_max(n)=min(real(x(m))+dist2anchor(n,m)+3*sigma,real(x_max(n)))...
           +j*min(imag(x(m))+dist2anchor(n,m)+3*sigma,imag(x_max(n)));
           % Do not cross borders of deployment area:
           if real(x_min(n))<min(area) x_min(n)=min(area)+j*imag(x_min(n)); end
           if imag(x_min(n))<min(area) x_min(n)=real(x_min(n))+j*min(area); end
           if real(x_max(n))>max(area) x_max(n)=max(area)+j*imag(x_max(n)); end
           if imag(x_max(n))>max(area) x_max(n)=real(x_max(n))+j*max(area); end
        else % if no 1 or 2 hop neighbors
            kk=kk+1;
            if kk==Na
               x_min(n)=min(area)+j*min(area); ; %complete area
               x_max(n)=max(area)+j*max(area); ; %complete area
            end
        end
    end
end

% Seed random variables (not necessary):
randn('state',0);
rand('twister',5489);

% Init. of beliefs:
for k=Na+1:Nn
	% Definition of variables:
	M_part(:,:,k)=zeros(Np,Niter);	% belief of target node k (particles) 
	M_wei(:,:,k)=zeros(Np,Niter);	% belief of target node k (weights) 
    M_bw(:,k)=zeros(Niter,1);      	% belief of target node k (bandwidth)
	% Initial state:
	M_part(:,1,k)=(real(x_max(k)-x_min(k))*rand(Np,1)+real(x_min(k)))...
    +j*(imag(x_max(k)-x_min(k))*rand(Np,1)+imag(x_min(k))); % initial particles (within b. box)
   	M_wei(:,1,k)=ones(Np,1)/Np;     % initial weights
    M_bw(1,k)=0.1+j*0.1;            	% initial bandwidth (bwx+j*bwy) (any value for first iter)
    % Initial belief:
    z(:,:,1,k)=ones(length(area),length(area)); % kde on grid (initially uniform)
end

% Init. of messages:
for m=Na+1:Nn
	for n=Na+1:Nn
        if dist(m,n)>0
            % Definition of variables:
        	msg_part(:,:,m,n)=zeros(Np, Niter);     % message from node m to node n (particles)
        	msg_wei(:,:,m,n)=zeros(Np, Niter);      % message from node m to node n (weights)
            msg_bw(:,m,n)=zeros(Niter,1);           % message from node m to node n (bandwidth)
            % Initial state:
        	msg_part(:,1,m,n)=(real(x_max(n)-x_min(n))*rand(Np,1)+real(x_min(n)))...
            +j*(imag(x_max(n)-x_min(n))*rand(Np,1)+imag(x_min(n))); % initial particles (within b. box)
            msg_wei(:,1,m,n)=ones(Np,1)/Np;         % initial weights            
            msg_bw(1,m,n)=0.1+j*0.1;                % initial bandwidth (bwx+j*bwy) (any value for first iter)
        end
	end
end

