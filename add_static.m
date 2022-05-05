function [Gs,unique_sta_num] = add_static( G, data )
disp('station static')
unique_sta_num = unique(data.ray.sta_num);
nstas = length( unique_sta_num );

sta_static = spalloc( size(G,1), nstas, size(G,1) ); 

for ista = 1:nstas
	sta_static(:,ista) = -( data.ray.sta_num == unique_sta_num(ista) );
end

Gs = [G sta_static];