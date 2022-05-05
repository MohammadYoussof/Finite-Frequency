%19 Jan 2012% M.Y.
function make_DS_map_synth_invention( mf, m, par, stax, stay)
% ds_map = 
disp(' CHECKERBOAD RESOLUTION TEST ya Mohammad :) ')

load data
stn = data.stn;
bbox = [min([stn.lon,stn.lat]); max([stn.lon,stn.lat])];
V = shaperead( 'unep_reported_malaria_cases_-_total_number_per_100000_population_world_1990-2003.shp', 'BoundingBox', bbox);
V1 =shaperead('ANG-15_boundaries.shp','BoundingBox', bbox);

Vcbf =shaperead('lips2011Polygons.shp','BoundingBox', bbox);
V3000 =shaperead('ne_10m_bathymetry_H_3000.shp','BoundingBox', bbox);

for i = 1:length(V3000)
	[x2,y2] = project_xy( par, V3000(i).Y, V3000(i).X );
	plot(x2,y2,'.','LineWidth',7,...
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

for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1, 'k' );
end
hold on
h = plot( stn.x, stn.y, 'ko' );
set( h, 'markerfacecolor', 'r' );

maxx = max( stn.x );
maxy = max( stn.y );
minx = min( stn.x );
miny = min( stn.y );
dx = (maxx - minx) * .1;
dy = (maxy - miny) * .1;
axis( [minx-dx maxx+dx miny-dy maxy+dy] );

cmap=colormap(redblue);
cmap=colormap(flipud(colormap));
colormap(cmap);

modx=par.modx; mody=par.mody; modz=par.modz; nx=length(modx); ny=length(mody); nz=length(modz);
plotx = modx(2):10:modx(end-1); ploty = mody(end-1):10:mody(2); plotz = modz(1:end-1);

% set limits for plot axes
min_x = min(plotx)+60;
max_x = max(plotx)-30;
min_y = min(ploty)+30;
max_y = max(ploty)-60;

for pz = 1:length(plotz)
    %clear m_plot
    m_plot1 = zeros(length(ploty),length(plotx));
    m_plot2 = zeros(length(ploty),length(plotx));
    for py = 1:length(ploty)
        for px = 1:length(plotx)
            pl_pt = [ ploty(py), plotx(px), plotz(pz) ];
            % find model node, upper-left-top
            my_diff = abs(mody - ploty(py)); srt_mdfy=sort(my_diff); mi = find(my_diff==srt_mdfy(1),1);
            mx_diff = abs(modx - plotx(px)); srt_mdfx=sort(mx_diff); mj = find(mx_diff==srt_mdfx(1),1);
            mz_diff = abs(modz - plotz(pz)); srt_mdfz=sort(mz_diff); mk = find(mz_diff==srt_mdfz(1),1);
            dx=0.5*(abs(modx(mj+1)-modx(mj))+abs(modx(mj-1)-modx(mj))); 
            dy=0.5*(abs(mody(mi)-mody(mi+1))+abs(mody(mi-1)-mody(mi)));     dh=sqrt(dx^2+dy^2);
            if mk>1
                dz=0.5*((modz(mk+1)-modz(mk))+(modz(mk)-modz(mk-1)));  
            elseif mk==1
                dz=modz(mk+1)-modz(mk);
            end
            dist = zeros(1,9);   val1 = zeros(1,9); val2 = zeros(1,9);
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
                        val1(jj+(ii-1)*3) = mf(mod_indx);
                        val2(jj+(ii-1)*3) = m(mod_indx);
                    %end
                end
            end
                wt = dist./sum(dist); 
                val1 = wt.*val1;
                val2 = wt.*val2;
                plot_val1 = -100*(sum(val1)/sum(wt));
                plot_val2 = -100*(sum(val2)/sum(wt));
                m_plot1(py,px) = plot_val1;
                m_plot2(py,px) = plot_val2;
        end
    end
    
                    % plot 1, output
    subplot(1,2,1);
    %m%  hold off
    
    
    hold on
    imagesc(plotx,-ploty,m_plot1)
    shading interp
     
    hold on
    plot( stax, stay, 'ko','MarkerFaceColor','k','MarkerSize', 3 ); 
     plot( x, y, 'k' );
     
    title( sprintf('Synthetic Input, Depth = %d\n',round(plotz(pz))) );
    colorbar
    caxis( [-0.6 0.6] );
    %view(0, 90)

    axis( [min_x max_x min_y max_y] );
    daspect([1 1 1])
    hold on
    
      
    xt=[-922 -700 -478 -256 -34 188 410 632 854];
    yt= [-772 -550 -328 -106 116 338 560 782];
 
    xt=[-922 -700 -478 -256 -34 188 410 632 854];
    yt= [-772 -550 -328 -106 116 338 560 782];
%     xt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
%     yt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
%     
%     xt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
%     yt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
    
    set(gca,'XTick',xt,'YTick',yt) 
    for gg=1:length(xt)
    [ytl,xtl] = project_xy( par, xt(gg)*ones(size(yt)), -yt, 'inverse' ); 
    xla(gg)=mean(xtl);
    end

    set(gca,'XTickLabel',roundn(xla,0.0),'YTickLabel',roundn(sort(ytl, 'descend'),0.0) )
    xlabel( 'Longitude ' );
    ylabel( 'Latitude ' );
        for i = 1:length(V3000)
	[x2,y2] = project_xy( par, V3000(i).Y, V3000(i).X );
	plot(x2,y2,'.','LineWidth',7,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
        end
hold on
    for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1, 'k' );
    end

    
    
                     % plot 2, output model %%%%
	subplot(1,2,2);
    
    
    load data
stn = data.stn;
bbox = [min([stn.lon,stn.lat]); max([stn.lon,stn.lat])];
V = shaperead( 'unep_reported_malaria_cases_-_total_number_per_100000_population_world_1990-2003.shp', 'BoundingBox', bbox);
for i = 1:length(V)
	[x,y] = project_xy( par, V(i).Y, V(i).X );
	plot( x, y, 'b' );
end
hold on
h = plot( stn.x, stn.y, 'ko' );
set( h, 'markerfacecolor', 'r' );

maxx = max( stn.x );
maxy = max( stn.y );
minx = min( stn.x );
miny = min( stn.y );
dx = (maxx - minx) * .1;
dy = (maxy - miny) * .1;
axis( [minx-dx maxx+dx miny-dy maxy+dy] );

    
    hold on
    
    %%hold off
    imagesc(plotx,-ploty,m_plot2)
    shading interp
    
     hold on
    plot( stax, stay, 'ko','MarkerFaceColor','k','MarkerSize', 3 ); 
     plot( x, y, 'k' );
    
    title( sprintf('Synthetic Result, Depth = %d\n',round(plotz(pz))) );
    colorbar
    caxis( [-0.6 0.6] );
    
    axis( [min_x max_x min_y max_y] );
    daspect([1 1 1])
    hold on
    
    %view(0, 90)
    %%%%plot( stax, -stay, 'ko','MarkerFaceColor','k','MarkerSize', 2 ); 
    
    
    
%     plot( Vx, -Vy, 'k' );  
%     xt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
%     yt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
     xt=[-922 -700 -478 -256 -34 188 410 632 854];
     yt= [-772 -550 -328 -106 116 338 560 782];

     xt=[-922 -700 -478 -256 -34 188 410 632 854];
     yt= [-772 -550 -328 -106 116 338 560 782];
%     xt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
%     yt=[-666 -555 -444 -333  -222 -111 0 111  222  333  444 555 666 ];
    
%     xt=[-776 -556 -337 -117 102 322 541 761]; yt= [-550 -350 -150 50 250 450 650];
%     xt=[-776 -556 -337 -117 102 322 541 761]; yt= [-542 -319 -96 125 349 571];
    set(gca,'XTick',xt,'YTick',yt) 
    for gg=1:length(xt)
    [ytl,xtl] = project_xy( par, -xt(gg)*ones(size(yt)), -yt, 'inverse' ); 
    xla(gg)=mean(xtl);
    end

    set(gca,'XTickLabel',roundn(xla,0.0),'YTickLabel',roundn(sort(ytl, 'descend'),0.0) )
    xlabel( 'Longitude ' );
    ylabel( 'Latitude ' );
    hold on
    for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1, 'k' );
    end
    hold on
    for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1, 'k' );
    end
hold on
    for i = 1:length(V3000)
	[x2,y2] = project_xy( par, V3000(i).Y, V3000(i).X );
	plot(x2,y2,'.','LineWidth',7,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
end
%hold on

% for i = 1:length(Vcbf)
% 	[x3,y3] = project_xy( par, Vcbf(i).Y, Vcbf(i).X );
% 	plot(x3,y3,'--rs','LineWidth',1,...
%                        'MarkerEdgeColor','k',...
%                        'MarkerFaceColor','y',...
%                        'MarkerSize',1);
%  
% end
    
%     choice = input( 'Print?', 's' );
% 	if (~ isempty(choice) )
            %print(['figures/Pwa_blocking1_' num2str(dh) '_' num2str(plotz(pz)) '.pdf'],'-dpdf','-r400');
% 	end
    pause(1)
    disp('Making Tomo maps is DONE ! :) YA MOHAMMAD') 
end
