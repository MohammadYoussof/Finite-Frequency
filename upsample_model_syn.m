function ds_map = upsample_model_syn( mf, m, par, type)
cmap=colormap(redblue);
cmap=colormap(flipud(colormap));
colormap(cmap);

%upsample_model(vm(:,inv),par,stax,stay);
modx=par.modx; mody=par.mody; modz=par.modz; nx=length(modx); ny=length(mody); nz=length(modz);
plotx = modx(2):10:modx(end-1); ploty = mody(end-1):10:mody(2); plotz = modz(2:end-1);

plot_z = modz(2:end-1); %modz(2):10:modz(end-1);
plot_z = setdiff(plot_z,modz);

%{
plot_pt.x   = zeros(length(plotx)*length(ploty)*(length(plot_z) + length(plotz) ),1 );
plot_pt.y   = zeros(size(plot_pt.x));
plot_pt.z   = zeros(size(plot_pt.x));
plot_pt.lat = zeros(size(plot_pt.x));
plot_pt.lon = zeros(size(plot_pt.x));
plot_pt.val = zeros(size(plot_pt.x));
%}
for ii = 1:(length(plot_z)+length(plotz))
    plot_pt_syn(ii).x   = zeros(length(plotx)*length(ploty),1);
    plot_pt_syn(ii).y   = zeros(size(plot_pt_syn(1).x));
    plot_pt_syn(ii).z   = zeros(size(plot_pt_syn(1).x));
    plot_pt_syn(ii).lat = zeros(size(plot_pt_syn(1).x));
    plot_pt_syn(ii).lon = zeros(size(plot_pt_syn(1).x));
    plot_pt_syn(ii).val = zeros(size(plot_pt_syn(1).x));
end

tic


clf
%n=1;
q=1;
for pz = 1:length(plotz)
    n=1;
    for py = 1:length(ploty)
        for px = 1:length(plotx)
            pl_pt_syn = [ ploty(py), plotx(px), plotz(pz) ];
            % find model node, upper-left-top
            my_diff = abs(mody - ploty(py)); srt_mdfy=sort(my_diff); mi = find(my_diff==srt_mdfy(1),1);
            mx_diff = abs(modx - plotx(px)); srt_mdfx=sort(mx_diff); mj = find(mx_diff==srt_mdfx(1),1);
            mz_diff = abs(modz - plotz(pz)); srt_mdfz=sort(mz_diff); mk = find(mz_diff==srt_mdfz(1),1);
            dx=0.5*(abs(modx(mj+1)-modx(mj))+abs(modx(mj-1)-modx(mj))); 
            dy=0.5*(abs(mody(mi)-mody(mi+1))+abs(mody(mi-1)-mody(mi)));     dh=1.15*sqrt(dx^2+dy^2);
            dist = zeros(1,9);   val2 = zeros(1,9);
            for ii = 1:3
                for jj = 1:3
                        mod_ptx = modx(mj+jj-2);
                        mod_pty = mody(mi+ii-2);
                        mod_ptz = modz(mk);
                        mod_pt = [ mod_pty mod_ptx mod_ptz ];
                        dist(jj+(ii-1)*3) = dh - sqrt( (mod_pt(1)-pl_pt_syn(1))^2 + (mod_pt(2)-pl_pt_syn(2))^2 );
                        if dist(jj+(ii-1)*3)<0
                            dist(jj+(ii-1)*3) = 0;
                        end
                        mi_ind = mi+ii-2; mj_ind = mj+jj-2; mk_ind = mk;
                        mod_indx = fast_sub2ind3d( [ny nx nz], mi_ind, mj_ind, mk_ind );
                        val2(jj+(ii-1)*3) = m(mod_indx); 
                end
            end
                wt = dist./sum(dist); 
                val2 = wt.*val2;
                
                plot_val2 = -100*(sum(val2)/sum(wt));
                plot_pt_syn(q).x(n) = pl_pt_syn(2);
                plot_pt_syn(q).y(n) = pl_pt_syn(1);
                plot_pt_syn(q).z(n) = pl_pt_syn(3);
                plot_pt_syn(q).val(n) = plot_val2;
                n=n+1;
        end
    end
    q=q+1;
end

for iz = 1:length(plot_z)
    mz_up = find(modz== max(modz(modz<plot_z(iz))));
    mz_dn = mz_up+1;
    up_dz = abs(plot_z(iz)-modz(mz_up));
    dn_dz = abs(modz(mz_dn)-plot_z(iz));
    dz    = modz(mz_dn)-modz(mz_up);
    up_wt = 1-(up_dz/dz);
    dn_wt = 1-(dn_dz/dz);
    plot_pt_syn(q).x   = plot_pt_syn(mz_up).x;
    plot_pt_syn(q).y   = plot_pt_syn(mz_up).y;
    plot_pt_syn(q).z   = plot_z(iz)*ones(size(plot_pt_syn(1).x));
    plot_pt_syn(q).val = (up_wt).*plot_pt_syn(mz_up).val + (dn_wt).*plot_pt_syn(mz_dn).val;
    q=q+1;
end

plot_pts_syn.x   = [];
plot_pts_syn.y   = [];
plot_pts_syn.z   = [];
plot_pts_syn.val = [];
for ii = 1:length(plot_z)+length(plotz)
    plot_pts_syn.x   = [ plot_pts_syn.x   ; plot_pt_syn(ii).x ];
    plot_pts_syn.y   = [ plot_pts_syn.y   ; plot_pt_syn(ii).y ];
    plot_pts_syn.z   = [ plot_pts_syn.z   ; plot_pt_syn(ii).z ];
    plot_pts_syn.val = [ plot_pts_syn.val ; plot_pt_syn(ii).val ];
end  
plot_pts_syn.lat = zeros(size(plot_pts_syn.x));
plot_pts_syn.lon = zeros(size(plot_pts_syn.x));

n=1:length(plot_pts_syn.x);
%for ii = 1:length(plot_pts.x)
    [plot_pts_syn.lat(n), plot_pts_syn.lon(n)] = project_xy(par, plot_pts_syn.x(n), plot_pts_syn.y(n), 'inverse');
%    n=n+1;
%end
if type==1
    save figures/plot_pts_syn plot_pts_syn
elseif type==2
    plot_pts_si = plot_pts_syn;
    save figures/plot_pts_si plot_pts_si
elseif type==3
    plot_pts_so = plot_pts_syn;
    save figures/plot_pts_so plot_pts_so
end
toc

