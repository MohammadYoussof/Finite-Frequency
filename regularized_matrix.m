function [GR,dR] = regularized_matrix( G, data, par, damp_param, ...
	smooth_param, LR )
 
d = data.ray.d;

R = [];
d_used = d;

LR = LR*smooth_param;

if ( par.use_crust_corr == true ) %use_calculated_static )
	d_used = d_used - data.ray.corr1;
end

[G_used,static_sta_num] = add_static( G, data );
[G_used, unique_orid]   = add_estatic( G_used, data );
num = size(G_used,2) - size(LR,2);
R = [R; [LR spalloc(size(LR,1),num,0)]];

par.damp_adjust = false;
if ( par.damp_adjust == true )
    G_damp = adjust_damp(G_used, data,par)*damp_param;
    R = [R; G_damp ]; 
else
    G_damp = reg_damp( data,par)*damp_param;
    R = [R; G_damp ]; 
end   

dR = [d_used; zeros(size(R,1),1)];
GR = [G_used; R];
