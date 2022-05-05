clear all
close all
clc
%reading and projecting map boundaries
load Namibia_P_Dataset_v03_W_normalVELOCITY_Sept_Oct_28
load data
for i=1:length(G)
m = G(i,:);
stn = data.stn;
%bbox = [min([stn.lon,stn.lat]); max([stn.lon,stn.lat])];
bbox=[[-87,-46];[120,8]];
V = shaperead( 'unep_reported_malaria_cases_-_total_number_per_100000_population_world_1990-2003.shp', 'BoundingBox', bbox);
V1 =shaperead('ANG-15_boundaries.shp','BoundingBox', bbox);
Vcbf =shaperead('lips2011Polygons.shp','BoundingBox', bbox);
V3000 =shaperead('ne_10m_bathymetry_H_3000.shp','BoundingBox', bbox);
for i = 1:length(V3000)
	[x2,y2] = project_xy( par, V3000(i).Y, V3000(i).X );
	plot(x2,y2,'.','LineWidth',10,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
end
hold on

for i = 1:length(Vcbf)
	[x3,y3] = project_xy( par, Vcbf(i).Y, Vcbf(i).X );
	plot(x3,y3,'--rs','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','y',...
                       'MarkerSize',1);
 
end

hold on
for i = 1:length(V)
	[x,y] = project_xy( par, V(i).Y, V(i).X );
	plot( x, y, 'k' );
end
hold on






h = plot( stn.x, stn.y, 'ko' );
set( h, 'markerfacecolor', 'b' );

maxx = max( stn.x );
maxy = max( stn.y );
minx = min( stn.x );
miny = min( stn.y );
dx = (maxx - minx) * .1;
dy = (maxy - miny) * .1;
axis( [minx-dx maxx+dx miny-dy maxy+dy] );

cmap=colormap(autumn);
cmap=colormap(flipud(colormap));
colormap(cmap);


modx=par.modx; mody=par.mody; modz=par.modz; nx=length(modx); ny=length(mody); nz=length(modz);
plotx = modx(2):10:modx(end-1); ploty = mody(end-1):10:mody(2); plotz = modz(1:end-1);
%plotx, ploty are modx and mody resampled to every 12 km and dropping the
%first and last nodes.

min_x = min(plotx)-300;
max_x = max(plotx)+300;
min_y = min(ploty)-300; 
max_y = max(ploty)+300;
%for i=10:10:100

clf
for pz = 1:length(plotz)
    %pz = 2:length(plotz)
    clear m_plot
    for py = 1:length(ploty)
        for px = 1:length(plotx)
            pl_pt = [ ploty(py), plotx(px), plotz(pz) ];
            % find model node, upper-left-top
            my_diff = abs(mody - ploty(py)); srt_mdfy=sort(my_diff); mi = find(my_diff==srt_mdfy(1),1);
            mx_diff = abs(modx - plotx(px)); srt_mdfx=sort(mx_diff); mj = find(mx_diff==srt_mdfx(1),1);
            mz_diff = abs(modz - plotz(pz)); srt_mdfz=sort(mz_diff); mk = find(mz_diff==srt_mdfz(1),1);
            dx=0.5*(abs(modx(mj+1)-modx(mj))+abs(modx(mj-1)-modx(mj))); 
            dy=0.5*(abs(mody(mi)-mody(mi+1))+abs(mody(mi-1)-mody(mi)));     dh= 1.15*sqrt(dx^2+dy^2);
            if mk>1
                dz=0.5*((modz(mk+1)-modz(mk))+(modz(mk)-modz(mk-1)));  
            elseif mk==1
                dz=modz(mk+1)-modz(mk);
            end
            dist = zeros(1,9);   val = zeros(1,9);
            for ii = 1:3
                for jj = 1:3
                    %for kk = 1:2
                        mod_ptx = modx(mj+jj-2);
                        mod_pty = mody(mi+ii-2);
                        mod_ptz = modz(mk);
                        mod_pt = [ mod_pty mod_ptx mod_ptz ];
                        dist(jj+(ii-1)*3) = dh - sqrt( (mod_pt(1)-pl_pt(1))^2 + (mod_pt(2)-pl_pt(2))^2 );
                        if dist(jj+(ii-1)*3)<0
                            dist(jj+(ii-1)*3) = 0;
                        end
                        mi_ind = mi+ii-2; mj_ind = mj+jj-2; mk_ind = mk;
                        mod_indx = fast_sub2ind3d( [ny nx nz], mi_ind, mj_ind, mk_ind );
                        val(jj+(ii-1)*3) = m(mod_indx); 
                    %end
                end
            end
                wt = dist./sum(dist); 
                val = wt.*val;
                plot_val = 100*(sum(val)/sum(wt));
                m_plot(py,px) = plot_val;
        end
    end
    %M% hold off
    
    imagesc(plotx,-ploty,m_plot)
    shading interp
    hold on
    plot( stax, -stay, 'ko','MarkerFaceColor','k','MarkerSize', 3 ); 
     plot( x, -y, 'k' );
    title( sprintf('Sensitivity Kernels , Depth Level= %f\n',plotz(pz)) );
    colorbar
    caxis( [-4.8 4.8] );
    
    %keyboard
    
    axis( [min_x max_x -max_y -min_y] ); %the wierdness with min_y, max_y is because we plot plotx, -ploty. See line 92
    daspect([1 1 1])


    xt=[-922 -700 -478 -256 -34 188 410 632 854];
    yt= [-772 -550 -328 -106 116 338 560 782];
 
    xt=[-922 -700 -478 -256 -34 188 410 632 854];
    yt= [-772 -550 -328 -106 116 338 560 782];
    
    set(gca,'XTick',xt,'YTick',yt) 
    for gg=1:length(xt)
    [ytl,xtl] = project_xy( par, xt(gg)*ones(size(yt)), -yt, 'inverse' ); 
    xla(gg)=mean(xtl);
    end

   % set(gca,'XTickLabel',roundn(xla,0),'YTickLabel',roundn(sort(ytl, 'descend'),0) )
   set(gca,'XTickLabel',roundn(xla,0.0),'YTickLabel',roundn(sort(ytl, 'descend'),0.0) )
   %keyboard
    xlabel( 'Longitude ' );
    ylabel( 'Latitude ' );
    drawnow
    %pause(1)
end
end