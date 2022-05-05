% ----------- resolution test ---------------------------------------------
if synth_test == true
    % ---- make a checkerboard pattern ------------
    mf = zeros( ny, nx, nz ); 
    [modx,mody,modz] = meshgrid( par.modx, par.mody, par.modz );

   amp = 0.015; 
    width = 2;
    rev_dep = 1;
    for iy = 1:floor(ny/width)
      indxy = ((iy-1)*width+2):iy*width;
      for ix = 1:floor(nx/width)
        indxx = ((ix-1)*width+2):ix*width;
        for iz = 1:floor(nz/rev_dep)
          indxz = ((iz-1)*rev_dep+1):iz*rev_dep;
          ischecker = false;
          if ( iseven(iz) )
            if ( (isodd(iy) && isodd(ix)) || (iseven(iy) && iseven(ix)) )
              ischecker = true;
            end
          else
            if ( (isodd(iy) && iseven(ix)) || (iseven(iy) && isodd(ix)) )
              ischecker = true;
            end
          end
          if ( ischecker )
            mf(indxy,indxx,indxz) = (amp);
          else
            mf(indxy,indxx,indxz) = (-amp);
          end
        end
      end 
    end

    % ---- invert synthetic dta ------------
    df = G*mf(:);
    DR = speye(nmodel);
    dR = [df; zeros(size(LR,1),1)];
    GR = [G; DR*damp];
    dR = [dR; zeros(nmodel,1)];
    GR = [GR; LR*smooth];
    [m,flag,relres,iter,resvec] = lsqr( GR, dR, 1e-6, 100 );
    % ---- plot synthetic inversion ------------
    make_DS_map_synth(mf,m,par,stax,stay); 

end