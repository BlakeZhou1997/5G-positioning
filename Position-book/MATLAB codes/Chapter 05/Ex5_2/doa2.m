function ydoa = doa2(r)
% Return the azimuth and elevation angles of a vector
ydoa=zeros(2,1);

ydoa(1) = atan2(r(2),r(1));
ydoa(2) = atan2(r(3),norm(r(1:2)));
