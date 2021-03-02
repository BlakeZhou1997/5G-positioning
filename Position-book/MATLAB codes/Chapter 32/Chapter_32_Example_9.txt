%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Chapter: Positioning in LTE
% Example: 9
% Title: Shape conversion from ellipsoid point with uncertainty ellipse
%        to ellipsoid point with uncertainty circle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% Parameters and variable setup
%

codedEllipsoidPoint = zeros(3,1); % Setup of array 
encodedEllipsoidPointWithUncertaintyCircle = zeros(1,5); % Setup of array
majorAxis = 6378137.0; % [m] Major axis of the WGS 84 Earth model, reference 3GPP TS 23.032, rev. 6.0.0. - see the UE Positioning FS
minorAxis = 6356752.314; % [m] Minor axis of the WGS 84 Earth model, reference 3GPP TS 23.032, rev. 6.0.0. - see the UE Positioning FS

%
% Tables according to Wigren et.al. 2010.
%

beta = zeros(32,1); % Initialize beta tables
gamma = zeros(32,8); % Intialize gamma tables

% Beta table

beta(:,1) = [
    1.0000;
    1.0771;
    1.1602;
    1.2496;
    1.3460;
    1.4497;
    1.5615;
    1.6819;
    1.8116;
    1.9513;
    2.1017;
    2.2638;
    2.4384;
    2.6264;
    2.8289;
    3.0470;
    3.2819;
    3.5350;
    3.8075;
    4.1011;
    4.4173;
    4.7579;
    5.1248;
    5.5200;
    5.9456;
    6.4040;
    6.8978;
    7.4296;
    8.0025;
    8.6195;
    9.2841;
   10.0000
];

% First gamma column 99% confidence

gamma(:,1) = [
    3.0349;
    2.9337;
    2.8543;
    2.7921;
    2.7445;
    2.7091;
    2.6837;
    2.6642;
    2.6488;
    2.6361;
    2.6268;
    2.6194;
    2.6124;
    2.6066;
    2.6028;
    2.5986;
    2.5957;
    2.5929;
    2.5900;
    2.5887;
    2.5868;
    2.5853;
    2.5838;
    2.5824;
    2.5818;
    2.5809;
    2.5807;
    2.5801;
    2.5789;
    2.5792;
    2.5788;
    2.5778
];

% Second gamma column 97% confidence

gamma(:,2) = [
    2.6488;
    2.5585;
    2.4833;
    2.4207;
    2.3707;
    2.3312;
    2.3002;
    2.2768;
    2.2571;
    2.2428;
    2.2315;
    2.2214;
    2.2135;
    2.2075;
    2.2019;
    2.1976;
    2.1932;
    2.1897;
    2.1867;
    2.1852;
    2.1829;
    2.1806;
    2.1793;
    2.1779;
    2.1772;
    2.1761;
    2.1755;
    2.1744;
    2.1738;
    2.1735;
    2.1735;
    2.1728
];

% Third gamma column 95% confidence

gamma(:,3) = [
    2.4480;
    2.3647;
    2.2919;
    2.2311;
    2.1801;
    2.1393;
    2.1054;
    2.0790;
    2.0583;
    2.0417;
    2.0277;
    2.0170;
    2.0082;
    2.0010;
    1.9951;
    1.9902;
    1.9852;
    1.9818;
    1.9789;
    1.9763;
    1.9737;
    1.9722;
    1.9705;
    1.9685;
    1.9683;
    1.9665;
    1.9663;
    1.9656;
    1.9642;
    1.9641;
    1.9633;
    1.9627
];

% Fourth gamma column 90% confidence

gamma(:,4) = [
    2.1468;
    2.0714;
    2.0054;
    1.9481;
    1.8986;
    1.8550;
    1.8192;
    1.7881;
    1.7636;
    1.7439;
    1.7276;
    1.7142;
    1.7034;
    1.6939;
    1.6864;
    1.6808;
    1.6758;
    1.6711;
    1.6677;
    1.6644;
    1.6610;
    1.6595;
    1.6577;
    1.6555;
    1.6539;
    1.6528;
    1.6521;
    1.6508;
    1.6497;
    1.6490;
    1.6485;
    1.6482
];

% Fifth gamma column 75% confidence

gamma(:,5) = [
    1.6659;
    1.6061;
    1.5526;
    1.5036;
    1.4588;
    1.4177;
    1.3810;
    1.3484;
    1.3193;
    1.2933;
    1.2713;
    1.2528;
    1.2365;
    1.2230;
    1.2113;
    1.2019;
    1.1949;
    1.1879;
    1.1828;
    1.1775;
    1.1740;
    1.1710;
    1.1676;
    1.1657;
    1.1631;
    1.1619;
    1.1600;
    1.1584;
    1.1579;
    1.1567;
    1.1556;
    1.1547
];

% Sixth gamma column 50% confidence

gamma(:,6) = [
    1.1780;
    1.1353;
    1.0958;
    1.0579;
    1.0228;
    0.9901;
    0.9595;
    0.9307;
    0.9035;
    0.8787;
    0.8560;
    0.8341;
    0.8149;
    0.7972;
    0.7806;
    0.7662;
    0.7527;
    0.7420;
    0.7320;
    0.7236;
    0.7157;
    0.7102;
    0.7051;
    0.7002;
    0.6966;
    0.6933;
    0.6912;
    0.6892;
    0.6863;
    0.6856;
    0.6839;
    0.6824
];

% Seventh gamma column 25% confidence

gamma(:,7) = [
    0.7594;
    0.7314;
    0.7050;
    0.6801;
    0.6565;
    0.6341;
    0.6116;
    0.5910;
    0.5711;
    0.5528;
    0.5349;
    0.5182;
    0.5015;
    0.4869;
    0.4720;
    0.4578;
    0.4453;
    0.4323;
    0.4208;
    0.4107;
    0.4009;
    0.3914;
    0.3821;
    0.3750;
    0.3671;
    0.3613;
    0.3547;
    0.3501;
    0.3457;
    0.3413;
    0.3379;
    0.3357
];

% Eight gamma column 1% confidence

gamma(:,8) = [
    0.1428;
    0.1378;
    0.1320;
    0.1281;
    0.1233;
    0.1178;
    0.1140;
    0.1094;
    0.1062;
    0.1023;
    0.0986;
    0.0951;
    0.0919;
    0.0877;
    0.0849;
    0.0821;
    0.0794;
    0.0759;
    0.0734;
    0.0710;
    0.0687;
    0.0654;
    0.0632;
    0.0610;
    0.0588;
    0.0567;
    0.0546;
    0.0525;
    0.0514;
    0.0493;
    0.0473;
    0.0462
];


% Step 1 - decoding

codedShape = positionFormatDescription(1); % The 3GPP code for the particular shape
codedEllipsoidPoint(1,1) = positionFormatDescription(2); % Coded latitude sign of the ellipsoid point, spare fields according to TS 23.032 are disregarded
codedEllipsoidPoint(2,1) = positionFormatDescription(3); % Coded degrees of latitude of the ellipsoid point
codedEllipsoidPoint(3,1) = positionFormatDescription(4); % Coded longitude of the ellipsoid point 
codedUncertaintySemiMajor = positionFormatDescription(5); % Coded semi-major axis of the uncertainty ellipse
codedUncertaintySemiMinor = positionFormatDescription(6); % Coded semi-major axis of the uncertainty ellipse
codedOrientationSemiMajor = positionFormatDescription(7); % Coded orientation of the semi major axis
if (length(positionFormatDescription)==8)
    codedConfidence = positionFormatDescription(8); % Coded confidence, optionality test added!
end
if (exist('codedConfidence'))
    confidenceBeforeTransformation = min([codedConfidence/100 1.0]); % Transform from % to fraction and override any configured passed value
end
uncertaintySemiMajor = 10*(1.1^codedUncertaintySemiMajor-1);
uncertaintySemiMinor = 10*(1.1^codedUncertaintySemiMinor-1);
    
% Step 2 - Confidence scaling

%
% Change against the limiting values of the confidence
%

if ( (confidenceBeforeTransformation==0.0) & (confidenceAfterTransformation==0.0) ) % First limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation==0.0) & (confidenceAfterTransformation>0.0) & (confidenceAfterTransformation<1.0) ) % Second limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation==0.0) & (confidenceAfterTransformation==1.0) ) % Third limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation>0.0) & (confidenceBeforeTransformation<1.0) & (confidenceAfterTransformation==0.0) ) % Fourth limiting case
    newUncertaintySemiMajor = 0.0; % Generating minimum coded value of 0...
    newUncertaintySemiMinor = 0.0; % % Generating minimum coded value of 0...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation>0.0) & (confidenceBeforeTransformation<1.0) & (confidenceAfterTransformation>0.0) & (confidenceAfterTransformation<1.0) ) % Fifth case - normal one, transform to nominal confidence shape
    scaleFactor2Nominal = 1/sqrt(-2*log(1-confidenceBeforeTransformation)); % Scale factor to get covariance confidence
    newUncertaintySemiMajor = scaleFactor2Nominal*uncertaintySemiMajor;
    newUncertaintySemiMinor = scaleFactor2Nominal*uncertaintySemiMinor;
    newConfidence = confidenceAfterTransformation;
elseif ( (confidenceBeforeTransformation>0.0) & (confidenceBeforeTransformation<1.0) & (confidenceAfterTransformation==1.0) ) % Sixth limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation==1.0) & (confidenceAfterTransformation==0.0) ) % Seventh limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation==1.0) & (confidenceAfterTransformation>0.0) & (confidenceAfterTransformation<1.0) ) % Eight limiting case
    newUncertaintySemiMajor = 0.0; % Generating minimum coded value of 0...
    newUncertaintySemiMinor = 0.0; % % Generating minimum coded value of 0...
    newConfidence = 0.0; % No confidence
elseif ( (confidenceBeforeTransformation==1.0) & (confidenceAfterTransformation==1.0) ) % Ninth limiting case
    newUncertaintySemiMajor = 2000000.0; % Generating maximum coded value of 127...
    newUncertaintySemiMinor = 2000000.0; % % Generating maximum coded value of 127...
    newConfidence = 0.0; % No confidence
else
end

% Step 3 - double interpolation for computation of uncertainty
% code

if (newUncertaintySemiMinor>newUncertaintySemiMajor) % Additional safety net
    temp = newUncertaintySemiMinor;
    newUncertaintySemiMinor = newUncertaintySemiMajor;
    newUncertaintySemiMajor = temp;
end

% The two interpolating tables are then found

if ( newConfidence>0.99)
    tableIndexLow = 1;
    tableIndexHigh = 1; % Saturate
    confidenceLow = 0.99;
    confidenceHigh = 100; % Avoids interpolating problems
elseif ( (newConfidence>0.97) & (newConfidence<=0.99) )
    tableIndexLow = 2;
    tableIndexHigh = 1;
    confidenceLow = 0.97;
    confidenceHigh = 0.99;
elseif ( (newConfidence>0.95) & (newConfidence<=0.97) )
    tableIndexLow = 3;
    tableIndexHigh = 2;
    confidenceLow = 0.95;
    confidenceHigh = 0.96;
elseif ( (newConfidence>0.90) & (newConfidence<=0.95) )
    tableIndexLow = 4;
    tableIndexHigh = 3;
    confidenceLow = 0.90;
    confidenceHigh = 0.95;
elseif ( (newConfidence>0.75) & (newConfidence<=0.90) )
    tableIndexLow = 5;
    tableIndexHigh = 4;
    confidenceLow = 0.75;
    confidenceHigh = 0.90;
elseif ( (newConfidence>0.50) & (newConfidence<=0.75) )
    tableIndexLow = 6;
    tableIndexHigh = 5;
    confidenceLow = 0.50;
    confidenceHigh = 0.75;
elseif ( (newConfidence>0.25) & (newConfidence<=0.50) )
    tableIndexLow = 7;
    tableIndexHigh = 6;
    confidenceLow = 0.25;
    confidenceHigh = 0.50;
elseif ( (newConfidence>0.01) & (newConfidence<=0.25) )
    tableIndexLow = 8;
    tableIndexHigh = 7;
    confidenceLow = 0.01;
    confidenceHigh = 0.25;
elseif (newConfidence<=0.01)
    tableIndexLow = 8;
    tableIndexHigh = 8; % Saturate
    confidenceLow = 0.01;
    confidenceHigh = 100; % Avoid interpolatimng problems
else
end

% The transformation to uncertainty circle then begins

if ( (newUncertaintySemiMinor==0.0) & (newUncertaintySemiMajor==0.0) ) % May be the result of a R99 report
    uncertaintyRadius = 0.0; 
    newConfidence = 0.0;
    uncertaintyCode = min([127 ceil(log(uncertaintyRadius/10+1)/log(1.1))]);
elseif ( (newUncertaintySemiMinor==2000000.0) & (newUncertaintySemiMajor==2000000.0) ) % May be the result of a R99 report
    uncertaintyRadius = 10*(1.1^127-1); 
    newConfidence = 0;
    uncertaintyCode = min([127 ceil(log(uncertaintyRadius/10+1)/log(1.1))]);
elseif ( (newUncertaintySemiMinor==0.0) & (newUncertaintySemiMajor>0.0) ) % Degenerate, handled by beta truncation
    thisBetaIndex = 32;
    gammaHigh = gamma(thisBetaIndex,tableIndexHigh); 
    gammaLow = gamma(thisBetaIndex,tableIndexLow); 
    thisGamma = ((gammaHigh-gammaLow)/(confidenceHigh-confidenceLow))*(newConfidence-confidenceLow)+gammaLow;
    uncertaintyRadius = thisGamma*newUncertaintySemiMajor;
    uncertaintyCode = min([127 ceil(log(uncertaintyRadius/10+1)/log(1.1))]);
else % Main case
    thisBeta = newUncertaintySemiMajor/newUncertaintySemiMinor;
    
    % Then find the beta indices corresponding to the beta value - loop
    % search need to be performed
    
    thisBetaIndexHigh = 1; % Initialization of search
    for i=1:32
        if ( beta(thisBetaIndexHigh,1)<thisBeta )
            thisBetaIndexHigh = thisBetaIndexHigh+1;
        end
    end
    if ( thisBetaIndexHigh==1 )
        thisBetaIndexLow = 1; % Saturate
    elseif ( thisBetaIndexHigh==33 )
        thisBetaIndexLow = 32;
        thisBetaIndexHigh = 32; % Saturate
    else
        thisBetaIndexLow = thisBetaIndexHigh-1;
    end
    
    % Intra-table interpolation first gives
    
    if ( thisBetaIndexLow==thisBetaIndexLow )
        gammaHigh = gamma(thisBetaIndexHigh,tableIndexHigh);
        gammaLow = gamma(thisBetaIndexLow,tableIndexLow);
    else
        gammaHigh = ((gamma(thisBetaIndexHigh,tableIndexHigh)-gamma(thisBetaIndexLow,tableIndexHigh))/(beta(thisBetaIndexHigh,1)-beta(thisBetaIndexLow,1)))*(thisBeta-beta(thisBetaIndexLow,1))+gamma(thisBetaIndexLow,tableIndexHigh);
        gammaLow = ((gamma(thisBetaIndexHigh,tableIndexLow)-gamma(thisBetaIndexLow,tableIndexLow))/(beta(thisBetaIndexHigh,1)-beta(thisBetaIndexLow,1)))*(thisBeta-beta(thisBetaIndexLow,1))+gamma(thisBetaIndexLow,tableIndexLow);
    end
    
    % Inter-table interpolation then follows
    
    thisGamma = ((gammaHigh-gammaLow)/(confidenceHigh-confidenceLow))*(newConfidence-confidenceLow)+gammaLow;
    uncertaintyRadius = thisGamma*newUncertaintySemiMajor;
    uncertaintyCode = min([127 ceil(log(uncertaintyRadius/10+1)/log(1.1))]);
end

% Step 4 - Encoding with uncertainty information

newCodedShape = 1; % Coded shape for ellipsoid point with uncertainty circle
encodedEllipsoidPointWithUncertaintyCircle(1) = newCodedShape;
encodedEllipsoidPointWithUncertaintyCircle(2) = codedEllipsoidPoint(1,1);
encodedEllipsoidPointWithUncertaintyCircle(3) = codedEllipsoidPoint(2,1);
encodedEllipsoidPointWithUncertaintyCircle(4) = codedEllipsoidPoint(3,1);
encodedEllipsoidPointWithUncertaintyCircle(5) = uncertaintyCode;


