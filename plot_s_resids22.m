function plot_s_resids22(par, data, s_static)

%load tomo_cmap;
%colormap( tomo_cmap );

colormap('default');
% S = shaperead('usastatehi');
% indx = find( strcmp({S.Name},'California') );
% [calx,caly] = project_xy( par, S(indx).Y, S(indx).X );
% indx2 = find( strcmp({S.Name},'Nevada') );
% [nevx,nevy] = project_xy( par, S(indx2).Y, S(indx2).X );
% indx3 = find( strcmp({S.Name},'Oregon') );
% [orx,ory] = project_xy( par, S(indx3).Y, S(indx3).X );
% indx4 = find( strcmp({S.Name},'Arizona') );
% [azx,azy] = project_xy( par, S(indx4).Y, S(indx4).X );
% indx5 = find( strcmp({S.Name},'Washington') );
% [wax,way] = project_xy( par, S(indx5).Y, S(indx5).X );
% indx6 = find( strcmp({S.Name},'Idaho') );
% [idx,idy] = project_xy( par, S(indx6).Y, S(indx6).X );

r1=stn.seen;
r2=stn.saad;
r3=stn.rawah;

min_x = min(r1)-100;
max_x = max(r1)+100;
min_y = min(r2)-100;
max_y = max(r2)+100;


clf
scatter( data.stn.x, data.stn.y, 400, s_static, '.' ); 
hold on
colorbar
caxis( [-0.5  0.5] );
axis( [min_x max_x min_y max_y] );
daspect([1 1 1])

% plot( calx, caly, 'k' );
% plot( nevx, nevy, 'k' );
% plot( orx, ory, 'k' );
% plot( azx, azy, 'k' );
% plot( wax, way, 'k' );
% plot( idx, idy, 'k' );    

xt=[]; yt=[];
set(gca,'XTick',xt,'YTick',yt) 


%xlabel( 'Longitude ' );
%ylabel( 'Latitude ' );
title('Station residuals (s)');
drawnow
%choice = input( 'Print?', 's' );
%if (~ isempty(choice) )
       print(['figures/stn_resids_try.ps'],'-depsc','-r300');
%end

% 
% figure(99)
% contourf(data.stn.x, data.stn.y, s_static)
%     cmap=colormap(jet); colormap(flipud(cmap));
%     set(gca,'clim', [-0.5 0.5])
%     colorbar

pause(1)