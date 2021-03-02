% Chapter 20
% -------------------------------------------------------------------------
% Example 1: Evaluation of the GDOP
% -------------------------------------------------------------------------
% author: F. Dovis

clear all; close all; clc

% Transmitters' positions

TX1=[-100,50];
TX2=[0, 50];
TX3=[100,0];
P=[0,0];

% Euclidean distances of the transmitters with respect to P

r1=sqrt((TX1(1)-P(1))^2+(TX1(2)-P(2))^2);     %TX1
r2=sqrt((TX2(1)-P(1))^2+(TX2(2)-P(2))^2);     %TX2
r3=sqrt((TX3(1)-P(1))^2+(TX3(2)-P(2))^2);     %TX3

% 1. Synchronous case

% Geometrical matrix
disp('User receiver synchronized to TX time scale')
Hs=[(TX1(1)-P(1))/r1, (TX1(2)-P(2))/r1;
    (TX2(1)-P(1))/r2, (TX2(2)-P(2))/r2;
    (TX3(1)-P(1))/r3, (TX3(2)-P(2))/r3]

HsTHs = Hs'*Hs;

rank_s = rank(HsTHs) %checks if the H matrix is singular

if rank_s>0
    Gs = inv(HsTHs);
    GDOP_s = sqrt(trace(Gs))
else
    disp('Not enough independent measurements available')
end

% 2. Asynchronous case

% Geometrical matrix
disp('User receiver not synchronized to TX time scale')
Hns=[(TX1(1)-P(1))/r1, (TX1(2)-P(2))/r1, -1;
    (TX2(1)-P(1))/r2, (TX2(2)-P(2))/r2, -1;
    (TX3(1)-P(1))/r3, (TX3(2)-P(2))/r3, -1]

HnsTHns = Hns'*Hns;

rank_ns = rank(HnsTHns) %checks if the H matrix is singular

if rank_ns>0
    Gns = inv(HnsTHns);
    GDOP_ns = sqrt(trace(Gns))
else
    disp('Not enough independent measurements available')
end

