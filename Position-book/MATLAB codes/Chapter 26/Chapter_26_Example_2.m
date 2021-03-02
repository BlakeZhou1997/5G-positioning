%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2.2: Draw particles within bounding box and compute message particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple example with 4 anchors and 1 unknown
% (see Example 2.3 for general version)
clear
Np=1000; % number of particles
sigma=0.1; % standard deviation of the measured distance (Gaussian noise)
area=[0 10]; % deployment area [min max] x [min max]

% True positions of 4 anchors (semi-random: near edges)
x(1)=rand(1,1)+j*rand(1,1);
x(2)=(rand(1,1)+max(area)-1)+j*(rand(1,1)+max(area)-1);
x(3)=rand(1,1)+j*(rand(1,1)+max(area)-1);
x(4)=(rand(1,1)+max(area)-1)+j*rand(1,1);
% Unknown node (semi-random position: in center)
x(5)=(max(area)/3-min(area)/3)*(rand(1,1)+j*rand(1,1))+max(area)/3*(1+j); 

% Measured distance (true+noise):
for k=1:4
    dist(k)=abs(x(5)-x(k)); % true distance
    dist(k)=dist(k)+3*sigma; % upper bound of the measured distance (in 99.7% cases)
end

% Bounds of the box (min-max method):
x_min=max(real(x(1:4))-dist(1:4))+j*max(imag(x(1:4))-dist(1:4));
x_max=min(real(x(1:4))+dist(1:4))+j*min(imag(x(1:4))+dist(1:4));
% Don't allow out of deployment area:
if real(x_min)<min(area) x_min=min(area)+j*imag(x_min); end
if imag(x_min)<min(area) x_min=real(x_min)+j*min(area); end
if real(x_max)>max(area) x_max=max(area)+j*imag(x_max); end
if imag(x_max)>max(area) x_max=real(x_max)+j*max(area); end

% Draw particles uniformly within the box:
particles=(real(x_max-x_min)*rand(Np,1)+real(x_min))+j*(imag(x_max-x_min)*rand(Np,1)+imag(x_min));

% Compute particles of the message of some neighboring unknown node:
dist2neigh=5+sigma*randn(Np,1); % noisy distance to neighboring unknown node
particles_msg=particles+dist2neigh.*exp(j*2*pi*rand(Np,1)); % shift in random direction
% Note: The weights of the messages are uniform

% Plot network and particles within box:
figure, 
plot(particles,'.g'); hold % plot initial particles
plot(real(x(1:4)),imag(x(1:4)),'sr','MarkerFaceColor','r'); % plot anchor nodes
plot(real(x(5)),imag(x(5)),'ok','MarkerFaceColor','k'); % plot unknown node
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]);hold off
% Plot network and particles from the message:
figure, 
plot(particles_msg,'.b'); hold % plot particles of the messages
plot(real(x(1:4)),imag(x(1:4)),'sr','MarkerFaceColor','r'); % plot anchor nodes
plot(real(x(5)),imag(x(5)),'ok','MarkerFaceColor','k'); % plot unknown node
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]);hold off

