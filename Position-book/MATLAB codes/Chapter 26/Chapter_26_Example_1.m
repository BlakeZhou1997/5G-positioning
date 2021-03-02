%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2.1: Compute pairwise potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
sigma=0.1; % standard deviation of the measured distance
xa=rand(1,1)+j*rand(1,1); % 2D position of the anchor
dist=0.3; % distance between anchor and unknown node
epsilon=1e-6; % some very small number
[X Y]=meshgrid(-0.5:0.02:1.5,-0.5:0.02:1.5); % example of 2D grid
pos=X+j*Y; % possible positions of unknown node (on grid)

logy=0.5*(abs(pos-xa)-dist).^2/sigma^2; % log of the potential
min_logy=min(min(logy)); % compute minimum (2 times for 2D grid)
y=exp(-(logy-min_logy)); % compute pairwise potential
% y=y.*exp(-min_logy) is not necessary (just a constant)
y=y+epsilon; % avoid zero
y=y/sum(sum(y)); % normalize
% 3D plot (pdf of unknown node, given anchor node):
mesh(X,Y,y);

