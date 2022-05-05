
function dv_out = dv_out(m,par)

modx=par.modx; mody=par.mody; modz=par.modz; nx=length(modx); ny=length(mody); nz=length(modz);

fid = fopen('dv_out_2','wt');
for pz = 1:nz
    for px = 1:nx
        for py = 1:ny
            mod_indx = fast_sub2ind3d( [ny nx nz], py, px, pz );
            x   = modx(px);
            y   = mody(py);        
            z   = modz(pz);
            val = -100*m(mod_indx); 
            fprintf(fid,'%d         %d          %d         %f\n',x,y,z,val);
        end
    end
end
fclose(fid);
