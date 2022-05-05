function indx = fast_sub2ind3d( dim, i, j, k )

if ( any(i > dim(1) | j > dim(2) | k > dim(3)) )
	error( 'TOMOLAB:FAST_SUB2IND3D:bad_index', ...
		'dimension index mismatch' );
end

indx = i + dim(1).*(j-1) + dim(1).*dim(2).*(k-1);
