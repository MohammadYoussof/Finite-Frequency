
%plotting resids map! %%%
function plot_s_resids_false(par, data, s_residuals)

load tomo_cmap;
colormap( tomo_cmap );

colormap('default');
%%%%pwave%%%%
stn = data.stn;
bbox = [min([stn.lon,stn.lat]); max([stn.lon,stn.lat])];
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
	plot( x, y, 'b' );
end
hold on
h = plot( stn.x, stn.y, 'ko' );
set( h, 'markerfacecolor', 'r' );

min_x = min(data.stn.stax)-100;
max_x = max(data.stn.stax)+100;
min_y = min(data.stn.stay)-100;
max_y = max(data.stn.stay)+100;


%%%%clf
scatter( data.stn.x, data.stn.y, 3000, s_residuals, '.' ); 
hold on

colorbar
caxis( [min(s_residuals)  max(s_residuals)] );
axis( [min_x max_x min_y max_y] );
daspect([1 1 1])

 %%%plot( x, y, 'k' );

xt=[]; yt=[];
set(gca,'XTick',xt,'YTick',yt) 


xlabel( 'Longitude ' );
ylabel( 'Latitude ' );
title('Station residuals (s)');
hold on

    for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, -y1, 'k' );
    end
    hold on
    
    for i = 1:length(V3000)
	[x2,y2] = project_xy( par, V3000(i).Y, V3000(i).X );
	plot(x2,-y2,'.','LineWidth',10,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
    end
hold on

    for i = 1:length(Vcbf)
	[x3,y3] = project_xy( par, Vcbf(i).Y, Vcbf(i).X );
	plot(x3,-y3,'--rs','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','y',...
                       'MarkerSize',1);
 
    end
    
box on 
grid on

drawnow
%choice = input( 'Print?', 's' );
%if (~ isempty(choice) )
       print(['figures_false_ref/Stns_Residuals_false.pdf'],'-dpdf','-r500');
%end


pause(1)