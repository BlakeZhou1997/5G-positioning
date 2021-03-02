function ydoa = doa(r)
% Return the azimuth and elevation angles of a vector
ydoa=zeros(1);

ydoa(1) = atan2(r(2),r(1));

