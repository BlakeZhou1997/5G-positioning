%Chapter 20
% -------------------------------------------------------------------------
% Example 2: Evaluation of the position
% -------------------------------------------------------------------------
% author:F. Dovis

clear all; close all; clc

%Transmitters' positions

TX1=[-100,50];
TX2=[0, 50];
TX3=[100,0];
P=[0,-10];

% Euclidean distances of the transmitters with respect to P

r1=sqrt((TX1(1)-P(1))^2+(TX1(2)-P(2))^2);     %TX1
r2=sqrt((TX2(1)-P(1))^2+(TX2(2)-P(2))^2);     %TX2
r3=sqrt((TX3(1)-P(1))^2+(TX3(2)-P(2))^2);     %TX3

% Measured pseudoranges (in meters)

rho1= 112.5;
rho2= 52;
rho3= 99;

% Differences between measurements and geometrical values

delta_rho=[r1-rho1; r2-rho2; r3-rho3]

%%%%%%%%%%%%%%%%%%%%%%%
% 1. Synchronous case %
%%%%%%%%%%%%%%%%%%%%%%%

% Geometrical matrix
disp('User receiver synchronized to TX time scale')
Hs=[(TX1(1)-P(1))/r1, (TX1(2)-P(2))/r1;
    (TX2(1)-P(1))/r2, (TX2(2)-P(2))/r2;
    (TX3(1)-P(1))/r3, (TX3(2)-P(2))/r3]

HsTHs = Hs'*Hs;
rank_s = rank(HsTHs) %checks if the H matrix is singular

if rank_s>0
    % Evaluation of the displacement wrt to the linearization point
    disp('Displacement wrt to the linearization point ')
    DeltaXs=inv(Hs'*Hs)*Hs'*delta_rho
    
    disp('Position ')
    Xs=P+DeltaXs'
else
    disp('Not enough independent measurements available')
end



%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Asynchronous case %
%%%%%%%%%%%%%%%%%%%%%%%%

% Geometrical matrix
disp('User receiver not synchronized to TX time scale')
Hns=[(TX1(1)-P(1))/r1, (TX1(2)-P(2))/r1, -1;
    (TX2(1)-P(1))/r2, (TX2(2)-P(2))/r2, -1;
    (TX3(1)-P(1))/r3, (TX3(2)-P(2))/r3, -1]

HnsTHns = Hns'*Hns;

rank_ns = rank(HnsTHns) %checks if the H matrix is singular

if rank_ns>0
    % Evaluation of the displacement wrt to the linearization point
    disp('Displacement wrt to the linearization point ')
    DeltaXns=inv(Hns'*Hns)*Hns'*delta_rho
    
    disp('Position ')
    Xns=P+DeltaXns(1:2)'
else
    disp('Not enough independent measurements available')
end



figure(1)
hold on
grid on
plot(TX1(1), TX1(2), 'o k','Linewidth',4)
plot(TX2(1), TX2(2), 'o b','Linewidth',4)
plot(TX3(1), TX3(2), 'o m','Linewidth',4)

plot(P(1), P(2), 's r','Linewidth',4)

plot(Xs(1),Xs(2), 'd y', 'Linewidth',4)
plot(Xns(1),Xns(2), 'd g', 'Linewidth',4)
legend('TX1', 'TX2', 'TX3', 'P', 'Xs' , 'Xns')

