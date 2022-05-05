function [raylat,raylon] = go_ray_trace( par, ray, stn, iray )

p = ray.p( iray );
baz = ray.baz( iray );
snum = ray.sta_num( iray );
slat = stn.lat( snum );
slon = stn.lon( snum );
v1 = par.vel.v1;
v2 = par.vel.v2;
z1 = par.vel.z1;
z2 = par.vel.z2;

[raylat,raylon] = ray_trace( v1, v2, z1, z2, p, slat, slon, ...
	baz );