function [Ges,unique_orid] = add_estatic( G, data )
disp('event static')
ray_orid = data.ray.orid;
unique_orid = unique( ray_orid );
norid = length( unique_orid );

evt_static = spalloc( size(G,1), norid, size(G,1) ); 

for iorid = 1:norid
	evt_static(:,iorid) = -( ray_orid == unique_orid(iorid) );
end

Ges = [G evt_static];