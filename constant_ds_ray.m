function rayxyz_out = constant_ds_ray( rayxyz, ds )

s = ray_length( rayxyz );
num = floor( s(end)/ds );
svec = [0:ds:num*ds]';

rayxyz_out = zeros( length(svec), 3 );

%{
rayxyz_out(:,1) = interp1( s, rayxyz(:,1), svec, 'linear*' );
rayxyz_out(:,2) = interp1( s, rayxyz(:,2), svec, 'linear*' );
rayxyz_out(:,3) = interp1( s, rayxyz(:,3), svec, 'linear*' );
%}

%{
rayxyz_out(:,1) = interp1q( s, rayxyz(:,1), svec );
rayxyz_out(:,2) = interp1q( s, rayxyz(:,2), svec );
rayxyz_out(:,3) = interp1q( s, rayxyz(:,3), svec );
%}

rayxyz_out(:,1) = interp1fast( s, rayxyz(:,1), svec );
rayxyz_out(:,2) = interp1fast( s, rayxyz(:,2), svec );
rayxyz_out(:,3) = interp1fast( s, rayxyz(:,3), svec );