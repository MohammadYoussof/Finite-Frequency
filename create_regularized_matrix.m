function [GR,dR] = create_regularized_matrix( G, data, par, damp, ...
	smooth, static_flag, LR )
disp('YES YEs Yes yes, it is now creating the regularized matrix')

use_inversion_static = ( strcmpi(static_flag,'invert') || ...
                        strcmpi(static_flag,'calc') );
use_calculated_static = ( strcmpi(static_flag,'calc') );

R = [];
d_used = data.ray.d;


LR = smooth*LR;


if ( use_calculated_static )
	d_used = d_used - data.ray.corr1;
end

if ( use_inversion_static )
	[G_used,static_sta_num] = add_static( G, data );
    [G_used, unique_orid]   = add_estatic( G_used, data );
    
    num = size(G_used,2) - size(LR,2);
    R = [R; [LR spalloc(size(LR,1),num,0)]];
else
	G_used = G;
    R = [R; LR];
end

G_damp = damp*speye(size(G_used,2));
R = [R; G_damp ];    

dR = [d_used; zeros(size(R,1),1)];
GR = [G_used; R];
