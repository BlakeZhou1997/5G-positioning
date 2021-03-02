function [stateAPosteriori,PAPostCovariance] =...
         KalmanRX(K, pdr, Upos, Uvel,...
         satPos, satVel, dt, ddt, AAA)

noSat = size(pdr,1); % number of available satellites
c = 299792458; % speed of light
f0 = 1.57542e9; % L1 carrier

% definition of the blocks of the measurement matrix
Range = zeros(noSat,1);
DeltaRange = zeros(noSat,1);
U = zeros(noSat,4);
for IndexU = 1 : noSat
    dX = (satPos(1,IndexU)*cos(AAA(IndexU)) +...
    satPos(2,IndexU)*sin(AAA(IndexU))) - Upos(1);
    dY = (satPos(2,IndexU)*cos(AAA(IndexU)) -...
    satPos(1,IndexU)*sin(AAA(IndexU))) - Upos(2);
    dZ 	= satPos(3,IndexU) - Upos(3);

    dXv     = satVel(1,IndexU) - Uvel(1);
    dYv     = satVel(2,IndexU) - Uvel(2);
    dZv     = satVel(3,IndexU) - Uvel(3);

    Range(IndexU) = sqrt(power(dX,2)+power(dY,2)+power(dZ,2));
    DeltaRange(IndexU) = ((dXv*dX)+(dYv*dY)+(dZv*dZ))/Range(IndexU);
    U(IndexU,1:3) = [dX dY dZ]/norm(Range(IndexU));
end
U(:,end) = -ones(noSat,1);

% definition of the measurements
PrCorr = pdr(:,2);
DeltaPrRaw = pdr(:,3)*c/f0;
DeltaPrCorr = -DeltaPrRaw + ddt;
zCurrentMeasurement = [PrCorr;DeltaPrCorr];

NominalMeasurement = [Range; DeltaRange];
KFMeasurement = zCurrentMeasurement - NominalMeasurement;

HObservationMatrix = [zeros(noSat*2,3),...
   kron([0 1; 1 0],-U(:,1:3)),zeros(noSat*2,8)];
HObservationMatrix(1:noSat,end-1) = -U(:,4);
HObservationMatrix(noSat+1:end-3,end) = -U(:,4);

[stateAPosteriori,PAPostCovariance] =...
    KalmanFilter(KFMeasurement,K.stateAPosteriori,...
    K.PAPostCovariance,HObservationMatrix,...
    K.RObservationNoiseCovariance,...
    K.PHIStateTransitMatrix,K.QStateNoiseCovariance);
