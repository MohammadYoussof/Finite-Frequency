function map_nodes = map_nodes( par, stax, stay )
clf
disp('Nodes Spacing Mapping')

% % %%%% JUST FOR NOW , I WILL MOVE IT LATER !!!%%%
% % 
 %=================== Input Rays ======================================%
    %fid = fopen('walpass2.file', 'r');
   fid = fopen('n3.dat', 'r');
   %%%%%% fid = fopen('offshore_sea.file', 'r');
    A = [textscan(fid, '%f %f %f %s %d %f %f %f %f %d'),1]; fclose(fid);
    ray.pd    = A{1};
    ray.p     = ray.pd .*(360/(2*pi*6371));
    ray.baz   = A{2};
    ray.d     = A{3};
    ray.sta   = A{4};
    ray.orid  = A{5};
    ray.cf    = A{6};
    ray.chan  = A{10};
    ray.lat   = A{7};
    ray.lon   = A{8};
    sta_ray   = ray.sta;
    ray.nrays = length(ray.p);

    %=================== Input Stations ==================================% 
    sta_names = unique( ray.sta );
    fid = fopen('new_crust_model.txt','r');
    B = [textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f'),1]; fclose(fid); 
    stns.sta = B{1};  stns.lat = B{2};  stns.lon = B{3};  stns.elv = B{4};
    stns.z   = B{5};  stns.n   = B{6};  stns.ne  = B{7};  stns.e   = B{8};
    stns.se  = B{9};  stns.s   = B{10}; stns.sw  = B{11}; stns.w   = B{12}; 
    stns.nw  = B{13}; 
    stn.lat = zeros(length(sta_names),1);
    stn.lon = zeros(length(sta_names),1);
    stn.elv = zeros(length(sta_names),1);
    stn.z   = zeros(length(sta_names),1);
    stn.n   = zeros(length(sta_names),1);
    stn.ne  = zeros(length(sta_names),1);
    stn.e   = zeros(length(sta_names),1);
    stn.se  = zeros(length(sta_names),1);
    stn.s   = zeros(length(sta_names),1);
    stn.sw  = zeros(length(sta_names),1);
    stn.w   = zeros(length(sta_names),1);
    stn.nw  = zeros(length(sta_names),1);
    stn.num = zeros(length(sta_names),1);
       %------ associate data -------------------%
    for ii = 1:length(sta_names)
        for jj = 1:length(stns.sta)
            if strcmp(sta_names(ii,:),stns.sta(jj,:))==1
                stn.lat(ii) = stns.lat(jj);
                stn.lon(ii) = stns.lon(jj);
                stn.elv(ii) = stns.elv(jj);
                stn.z(ii)   = stns.z(jj);
                stn.n(ii)   = stns.n(jj);
                stn.ne(ii)  = stns.ne(jj);
                stn.e(ii)   = stns.e(jj);
                stn.se(ii)  = stns.se(jj);
                stn.s(ii)   = stns.s(jj);
                stn.sw(ii)  = stns.sw(jj);
                stn.w(ii)   = stns.w(jj);
                stn.nw(ii)  = stns.nw(jj);
                stn.num(ii) = ii;
            end
        end
    end 

    stn.sta = sta_names;
    ray.sta_num = zeros(length(ray.d),1);
    for ii = 1:length(ray.d)
        for jj = 1:length(stn.sta)
            if strcmp(char(ray.sta(ii)),stn.sta(jj))==1
                ray.sta_num(ii) = jj;
            end
        end
    end 
    sta_names = char(sta_names);
    stn.sta   = sta_names;
    stn.nstas = length(sta_names);
    
    % data structure 
    data.stn = stn;
    data.ray = ray;
    


%%%%%%%%%%%%#
%data

disp 'here we go'
load data
stn = data.stn;
bbox = [min([stn.lon,stn.lat]); max([stn.lon,stn.lat])];
V = shaperead( 'unep_reported_malaria_cases_-_total_number_per_100000_population_world_1990-2003.shp', 'BoundingBox', bbox);
V1 =shaperead('ANG-15_boundaries.shp','BoundingBox', bbox);
for i = 1:length(V)
	[x,y] = project_xy( par, V(i).Y, V(i).X );
	plot( x, y );
end
%hold on

for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1 );
end
%stn.x=stax; stn.y=stay;
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

cmap=colormap(hot);
cmap=colormap(flipud(colormap));
colormap(cmap);

modx=par.modx; mody=par.mody; modz=par.modz; nx=length(modx); ny=length(mody); nz=length(modz);
plotx = modx(2):15:modx(end-1); ploty = mody(end-1):15:mody(2); plotz = modz(1:end-1);

 min_x = min(stax)-500;
 max_x = max(stax)+500;
 min_y = min(stay)-500;
 max_y = max(stay)+500;

if ( (max_x-min_x) >= (max_y-min_y) )
    p_min = min_x;
    p_max = max_x;
else
    p_min = min_y;
    p_max = max_y;
end
%%%%%------%%%%
%hold off
  pz=1;
    m_plot = zeros(length(ploty),length(plotx));
    for py = 1:length(ploty)
        for px = 1:length(plotx)
            my_diff = abs(mody - ploty(py)); srt_mdfy=sort(my_diff); mi = find(my_diff==srt_mdfy(1),1);
            mx_diff = abs(modx - plotx(px)); srt_mdfx=sort(mx_diff); mj = find(mx_diff==srt_mdfx(1),1);
            mz_diff = abs(modz - plotz(pz)); srt_mdfz=sort(mz_diff); 
            %mk = find(mz_diff==srt_mdfz(1),1);
            %dx=0.5*(abs(modx(mj+1)-modx(mj))+abs(modx(mj-1)-modx(mj))); 
            %dy=0.5*(abs(mody(mi)-mody(mi+1))+abs(mody(mi-1)-mody(mi)));     
            
            if ( mj>3 && mj<(length(modx)-2) && mi>3 && mi<(length(mody)-2) )
                            dhn = 0.25*( (modx(mj+1)-modx(mj))+(modx(mj)-modx(mj-1)) ...
                            +(mody(mi-1)-mody(mi))+(mody(mi)-mody(mi+1)) );
            elseif ( mj<=3 || mj>=(length(modx)-2) || mi<=3 || mi>=(length(mody)-2) )
                            dhn = modx(2)-modx(1);
            else
                            error('dh error');
            end             
            m_plot(py,px) = dhn;

        end
    end
    hold on
    imagesc(plotx,-ploty,m_plot)
    shading interp
    hold on
    plot( stax, stay, 'ko','MarkerFaceColor','k','MarkerSize', 3 ); 
     plot( x, y, 'k' );
    title( sprintf('Interpolated average node-spacing') );
    colorbar
    caxis( [30 45] );
    
    
    xt=p_min:150:p_max; yt=p_min:150:p_max;
    set(gca,'XTick',xt,'YTick',yt) 
    
    for gg=1:length(xt)
    [ytl,xtl] = project_xy( par, xt(gg)*ones(size(yt)), -yt, 'inverse' ); 
    xla(gg)=mean(xtl);
    end
%     
    set(gca,'XTickLabel',roundn(xla,0.0),'YTickLabel',roundn(sort(ytl, 'ascend'),0.0) )
    xlabel( 'Longitude' );
    ylabel( 'Latitude' );
	hold on
    for i = 1:length(V1)
	[x1,y1] = project_xy( par, V1(i).Y, V1(i).X );
	plot( x1, y1 );
    end
    drawnow

	%choice = input( 'Print?', 's' );
	%if (~ isempty(choice) )
          print('figures/Nodes_Mapping.pdf','-dpdf','-r500');
	%end
    pause(1)

