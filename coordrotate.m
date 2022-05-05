function [long,cross] = coordrotate(lat,lon) 

A=[17.00 -35.00];
B=[33.00 -18.00];

%A= [7.00 -18.60];
%B= [16.00 -18.50];
rotangle=atan((B(2)-A(2))/(B(1)-A(1)));
rotangle=abs(rotangle);
lat=lat-A(2);
lon=lon-A(1);

long=lon*cos(rotangle)-lat*sin(rotangle);
cross=lon*sin(rotangle)+lat*cos(rotangle);

% lon*cos(rotangle)
% lat*sin(rotangle)
% lon*sin(rotangle)
% lat*cos(rotangle)

%keyboard