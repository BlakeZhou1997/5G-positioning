%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Chapter: Positioning in LTE
% Example: 2
% Title: Computation of the QoS of an ellipse
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (positionFormat == 'ELLIPSOID_POINT_WITH_UNCERTAINTY_ELLIPSE')

    % Substep 1 - decoding of format
        
    % Unpacking of format
        
    codedShape = positionFormatDescription(1); % Coded shape
    codedEllipsoidPoint(1) = positionFormatDescription(2); % Latitude sign
    codedEllipsoidPoint(2) = positionFormatDescription(3); % Degrees of latitude
    codedEllipsoidPoint(3) = positionFormatDescription(4); % Longitude
    codedUncertaintySemiMajor = positionFormatDescription(5); % Semi-major axis
    codedUncertaintySemiMinor = positionFormatDescription(6); % Semi-minor axis
    codedOrientationSemimajor = positionFormatDescription(7); % Orientation od semi-major axis
    confidence = positionFormatDescription(8); % Confidence
        
    % Calculation of the relevant variables
        
    uncertaintySemiMajor = 10*(1.1^codedUncertaintySemiMajor-1);
    uncertaintySemiMinor = 10*(1.1^codedUncertaintySemiMinor-1);
        
    % Substep 2 - calculation of accuracy
        
    areaEllipsoidPointWithUncertaintyEllipse = pi*uncertaintySemiMajor*uncertaintySemiMinor; % Uncertainty area
    equivalentUncertaintyRadiusEllipsoidPointWithUncertaintyEllipse = sqrt(areaEllipsoidPointWithUncertaintyEllipse/pi);
    accuracyCodeEllipsoidPointWithUncertaintyEllipse = min([127 ceil(log(equivalentUncertaintyRadiusEllipsoidPointWithUncertaintyEllipse/10+1)/log(1.1))]);
       
    % Substep 4 - Reporting
        
    AccuracyCode = accuracyCodeEllipsoidPointWithUncertaintyEllipse;

end