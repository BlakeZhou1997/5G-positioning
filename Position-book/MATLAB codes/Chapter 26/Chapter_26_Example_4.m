%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3.1: Draw 6D particles from 3-node clique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple example with 4 anchors and one clique (3 unknowns)
clear
Np=1000; % number of particles
sigma=0.1; % standard deviation of the measured distance (Gaussian noise)
area=[0 10]; % deployment area [min max] x [min max]

% True positions of 4 anchors (semi-random: near edges)
x(1)=rand(1,1)+j*rand(1,1);
x(2)=(rand(1,1)+max(area)-1)+j*(rand(1,1)+max(area)-1);
x(3)=rand(1,1)+j*(rand(1,1)+max(area)-1);
x(4)=(rand(1,1)+max(area)-1)+j*rand(1,1);

% True positions of unknowns (random)
for k=5:7
    x(k)=(max(area)-min(area))*(rand(1,1)+j*rand(1,1));
end
% Uncomment following line if you want determinsitic placement:
%x(5)=1+j*1; x(6)=8+j*3; x(7)=6+j*9;

% Measured distance (true+noise):
% (fully connected network)
for m=1:7
    for n=5:7 % we don't need anchor-anchor distances       
        dist_true(m,n)=abs(x(m)-x(n)); % true distance
        dist(m,n)=dist_true(m,n)+3*sigma; % upper bound of the measured distance (in 99.7% cases)
    end
end

% Bounds of the boxes (min-max method):
for n=5:7
    x_min(n)=max(real(x(1:4))-dist(1:4,n)')+j*max(imag(x(1:4))-dist(1:4,n)');
    x_max(n)=min(real(x(1:4))+dist(1:4,n)')+j*min(imag(x(1:4))+dist(1:4,n)');
    % Don't let out of deployment area:
    if real(x_min(n))<min(area) x_min(n)=min(area)+j*imag(x_min(n)); end
    if imag(x_min(n))<min(area) x_min(n)=real(x_min(n))+j*min(area); end
    if real(x_max(n))>max(area) x_max(n)=max(area)+j*imag(x_max(n)); end
    if imag(x_max(n))>max(area) x_max(n)=real(x_max(n))+j*max(area); end
    % Draw particles uniformly within the box:
    particles_naive(:,n)=(real(x_max(n)-x_min(n))*rand(Np,1)+real(x_min(n)))...
        +j*(imag(x_max(n)-x_min(n))*rand(Np,1)+imag(x_min(n)));    
end

% Draw particles within the box and accept if satisfy distance restriction:
for p=1:Np % for each particle
    sim_max=1; 
    while (sim_max<100) % check max. 50 times
        % Draw candidate for the particle from each node:
        particles(p,5)=(real(x_max(5)-x_min(5))*rand(1,1)+real(x_min(5)))...
            +j*(imag(x_max(5)-x_min(5))*rand(1,1)+imag(x_min(5)));
        particles(p,6)=(real(x_max(6)-x_min(6))*rand(1,1)+real(x_min(6)))...
            +j*(imag(x_max(6)-x_min(6))*rand(1,1)+imag(x_min(6))); 
        particles(p,7)=(real(x_max(7)-x_min(7))*rand(1,1)+real(x_min(7)))...
            +j*(imag(x_max(7)-x_min(7))*rand(1,1)+imag(x_min(7)));
        % Distance must be mean (in our case, true) plus/minus 3*sigma:
        if not(abs(abs(particles(p,5)-particles(p,6))-dist_true(5,6))<3*sigma) ...
        	|| not(abs(abs(particles(p,5)-particles(p,7))-dist_true(5,7))<3*sigma)...
        	|| not(abs(abs(particles(p,6)-particles(p,7))-dist_true(6,7))<3*sigma)
         	  sim_max=sim_max+1; 
        else break; % accept particle
        end
    end
end


% Plot:
% Clique particles (naive):
figure, 
plot(particles_naive(:,5),'.g'); hold % plot particles of node 5
plot(particles_naive(:,6),'.m'); % plot particles of node 6
plot(particles_naive(:,7),'.b'); % plot particles of node 7
plot(real(x(1:4)),imag(x(1:4)),'sr','MarkerFaceColor','r'); % plot anchors
plot(real(x(5:7)),imag(x(5:7)),'ok','MarkerFaceColor','k'); % plot unknowns
plot([real(x(5)) real(x(6))],[imag(x(5)) imag(x(6))],'-k') % plot line
plot([real(x(5)) real(x(7))],[imag(x(5)) imag(x(7))],'-k') % plot line
plot([real(x(6)) real(x(7))],[imag(x(6)) imag(x(7))],'-k') % plot line
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]); hold off
% Clique particles (improved):
figure, 
plot(particles(:,5),'.g'); hold % plot particles of node 5
plot(particles(:,6),'.m'); % plot particles of node 6
plot(particles(:,7),'.b'); % plot particles of node 7
plot(real(x(1:4)),imag(x(1:4)),'sr','MarkerFaceColor','r'); % plot anchors
plot(real(x(5:7)),imag(x(5:7)),'ok','MarkerFaceColor','k'); % plot unknowns
plot([real(x(5)) real(x(6))],[imag(x(5)) imag(x(6))],'-k') % plot line
plot([real(x(5)) real(x(7))],[imag(x(5)) imag(x(7))],'-k') % plot line
plot([real(x(6)) real(x(7))],[imag(x(6)) imag(x(7))],'-k') % plot line
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]); hold off

