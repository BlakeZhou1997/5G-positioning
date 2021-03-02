%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Chapter: Positioning in LTE
% Example: 3
% Title: Coordinate transformations WGS84 <--> ET
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1 - Unpacking of format and decoding

majorAxis = 6378137.0; % WGS84 earth ellipsoid major axis
minorAxis = 6356752.314; % WGS84 earth ellipsoid minor axis
f = 1-minorAxis/majorAxis;
eSquared = 2*f-f*f; 

N_p=positionFormatDescription(2);

for i=1:N_p  % Cell Polygon corner unpacking
   CellPosition(3*i-2) = positionFormatDescription(3*i);
   CellPosition(3*i-1) = positionFormatDescription(3*i+1);
   CellPosition(3*i) = positionFormatDescription(3*i+2);
end

% Decoding of cell polygon information, to express corners in latitudes and
% longitudes

for i=1:N_p
    CellLatLong(1,i) = (1-2*CellPosition(3*i-2))*CellPosition(3*i-1)*(pi/2)/2^23; % i:th corner latitude 
    CellLatLong(2,i) = CellPosition(3*i)*(2*pi)/2^24; % i:th corner longitude
end

% Step 2 - Transformation to north-east earth tangential Cartesian system
        
for i=1:N_p % First transform to earth centered Cartesian coordinates
    N = majorAxis/sqrt(1-eSquared*(sin(CellLatLong(1,i)))^2); % Re-used variable in coordinate transformation 
    CellXYZ(1,i) = N*cos(CellLatLong(1,i))*cos(CellLatLong(2,i)); % Earth centered x-coordinate of i:th cell polygon corner
    CellXYZ(2,i) = N*cos(CellLatLong(1,i))*sin(CellLatLong(2,i)); % Earth centered y-coordinate of i:th cell polygon corner
    CellXYZ(3,i) = N*(minorAxis/majorAxis)^2*sin(CellLatLong(1,i)); % Earth centered z-coordinate of i:th cell poltgon corner, automatically on WGS 84 ellipsoid
end
for i=1:N_p % Then transform to earth tangential Cartesian north - east co-ordinates, rotation around the origin which is in the first corner of the cell polygon
    polygonXY(1,i) = -sin(CellLatLong(2,1))*(CellXYZ(1,i)-CellXYZ(1,1))+cos(CellLatLong(2,1))*(CellXYZ(2,i)-CellXYZ(2,1)); % Earth tangential x-coordinate of i:th polygon corner
    polygonXY(2,i) = -sin(CellLatLong(1,1))*cos(CellLatLong(2,1))*(CellXYZ(1,i)-CellXYZ(1,1))-sin(CellLatLong(1,1))*sin(CellLatLong(2,1))*(CellXYZ(2,i)-CellXYZ(2,1))+cos(CellLatLong(1,1))*(CellXYZ(3,i)-CellXYZ(3,1)); % Earth-tangential y-coordinate of the i:th polygon corner 
end
     
% Step 3 - Copy to new polygon for back transformation

polygonNewXY = polygonXY; % Normally something additional is done here, e.g. up/doen scaling computations.

% Step 4 - Coordinate transformation back to WGS84 latitude/longitude

%
% The first step is to go from the earth tangential system back to earth
% centered coordinates. Again, the first corner represents the origin. 
%

for i=1:N_p
    NewCellXYZ(1,i) = CellXYZ(1,1)-sin(CellLatLong(2,1))*polygonNewXY(1,i)-sin(CellLatLong(1,1))*cos(CellLatLong(2,1))*polygonNewXY(2,i);
    NewCellXYZ(2,i) = CellXYZ(2,1)+cos(CellLatLong(2,1))*polygonNewXY(1,i)-sin(CellLatLong(1,1))*sin(CellLatLong(2,1))*polygonNewXY(2,i);
    NewCellXYZ(3,i) = CellXYZ(3,1)+cos(CellLatLong(1,1))*polygonNewXY(2,i);
end

% Then transform to WGS 84 lat/long. Result is in radians

for i=1:N_p
    N = sqrt(NewCellXYZ(1,i)^2+NewCellXYZ(2,i)^2+(NewCellXYZ(3,i)*majorAxis^2/minorAxis^2)^2);
    NewCellLatSign(1,i) = sign(NewCellXYZ(3,i));
    NewCellLatLong(1,i) = NewCellLatSign(1,i)*acos(sqrt((NewCellXYZ(1,i)^2+NewCellXYZ(2,i)^2)/N^2));
    NewCellLatLong(2,i) = atan2(NewCellXYZ(2,i),NewCellXYZ(1,i));
end

% Step 5 - Polygon encoding. 

codedShape = 5; % 3GPP polygon code
encodedPolygon = [codedShape N_p]; % The first two IEs of the 3GPP polygon shape 

% Encode and add the corners to encodedPolygon one by one

for i=1:N_p
    if ( NewCellLatSign(1,i)>0.0 )
        thisCodedPolygonCornerPoint(1,1) = 0; % Code for North
    else
        thisCodedPolygonCornerPoint(1,1) = 1; % Code for South
    end
    if ( abs(NewCellLatLong(1,i))>= pi/2 ) % Additional safety net
        thisCodedPolygonCornerPoint(2,1) = 2^23-1;
    else
        thisCodedPolygonCornerPoint(2,1) = floor((2^24/pi)*abs(NewCellLatLong(1,i)));
    end
    if ( NewCellLatLong(2,i)>=pi )
        NewCellLatLong(2,i) = -pi;
    elseif ( NewCellLatLong(2,i)<-pi )
        NewCellLatLong(2,i) = -pi;
    end
    thisCodedPolygonCornerPoint(3,1) = floor((2^23/pi)*NewCellLatLong(2,i));
    
    encodedPolygon = [encodedPolygon thisCodedPolygonCornerPoint(1,1) thisCodedPolygonCornerPoint(2,1) thisCodedPolygonCornerPoint(3,1)]; % Dynamic addition of one polygon corner point
end
