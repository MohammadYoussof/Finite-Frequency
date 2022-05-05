function stds = look_delays( par, data, flag )
disp ('Checking DelayTimes...')
%load Namibia_P_Dataset_ref1.mat;
stn = data.stn;
ray = data.ray;


unique_sta_num = unique( ray.sta_num );
nsta = length( unique_sta_num );
stds = zeros( 1, nsta );

nplots_i = 2;
nplots_j = 2;
nplots = nplots_i * nplots_j;
for i = 1:nplots_i
	for j = 1:nplots_j
		
		indx = sub2ind( [nplots_i,nplots_j], i, j );
		h(i,j) = subplot(nplots_i,nplots_j, indx );
		
	end
end

% load tomo_cmap;
% colormap( tomo_cmap );
cmap = colormap(redblue);
cmap = flipud(cmap);
colormap(cmap);

for ista = 1:nplots:nsta
	
	for iplot = 1:nplots
		cla( h(iplot), 'reset' );
	end
	
	for iplot = 1:nplots
		
		sta_indx = ista+iplot-1;
        
	
		if ( sta_indx > nsta )
			break;
		end
		
		sta = stn.sta(sta_indx,:);
		
		indx = find(unique_sta_num(sta_indx) == ray.sta_num);
		baz = dfe2az( deg2rad(ray.baz(indx)) );
		inc = rad2deg(asin(ray.p(indx)*8));
		
		switch ( flag )
			
			case 'raw'
				
				d = ray.d(indx);
			
             case 'invented'
% 				load Namibia_P_Dataset_ref1.mat;
%                 dc1=dc;
%                 
 				d = dc1(indx);
            
			case 'corrected'
				
				if ( isfield( ray, 'corr1' ) )
					d = ray.d(indx) - ray.corr1(indx);
                  
				else
					error( 'TOMOLAB:LOOK_DELAYS:no_corrections', ...
						'no corrections' );
				end
				
			case 'corrections'
				
				d = ray.corr1(indx);
			
			case 'residuals'
				
				d = data.ray.resids(indx);
			
			otherwise
				
				error( 'PTOMO:look_delays:bad_flag', ...
					'flag must be ''corrected'', or ''corrections''' );
				
        end
          plot_delays( h(iplot), baz, inc, d,  sta, cmap, flag );
% % % % % % % % % % % % % % % % % % % % % % % % 		plot_delays( h(iplot), baz, inc, d,  sta, tomo_cmap, flag );
% % % %         save baz1 baz
% % % %         save inc1 inc
        %by mohammad /10/1/2012
        %m%save ( 'h(iplot)', 'baz', 'inc', 'd','sta', 'tomo_cmap')
        %m%save plot_delays
        
        
            
        
    end	
	pause
	drawnow
    %disp(['The plot is: ' num2str(ista)])
	 
end

%disp ('Checking DelayTimes'  num2str(ista))

save stds1 stds

function plot_delays( ax, baz, inc, d, sta, cmap, flag )

for gaz = 1:length(ax)
    
x = cos(baz) .* inc;
y = sin(baz) .* inc;
hold( ax, 'off' );
h = scatter( ax, x, y, 900, d, '.' );
hold( ax, 'on' );

circ = deg2rad(0:361);
x = cos( circ ) * 15;
y = sin( circ ) * 15;
plot( ax, x, y, 'k' );
%i uncommented this
text( x(1), y(1), '30^\circ' );

x = cos( circ ) * 30;
y = sin( circ ) * 30;
plot( ax, x, y, 'k' );
%i uncommented this
text( x(1), y(1), '60^\circ' );

x = cos( circ ) * 45;
y = sin( circ ) * 45;
plot( ax, x, y, 'k' );
%i uncommented this follwong line
text( x(1), y(1), '90^\circ' );

h = plot( ax, 0, 0, 'k+' );
set( h, 'markersize', 55 );

txt = sprintf( 'Station %s, %.2f to %.2f s, std = %.2f', ...
	sta, min(d), max(d), std(d) );
fid = fopen('delaytimes_model','a');
fprintf(fid , ' %s  %f  %f  %f  %f \n', sta, min(d), max(d), std(d), mean(d) )
fclose(fid);

title( ax, txt, 'fontsize', 10, 'fontweight', 'bold' );



axis( ax, 'equal' );
axis( ax, 'tight' );
grid( ax, 'on' );

switch (flag)
    case 'residuals'
        val = 1.00;
    otherwise
        val = 1.20;
end
%%%%%%MMMM%%%%% caxis( ax, [(-val) (val)] );
caxis( ax, [(-0.8) (0.8)] );
colorbar( 'peer', ax, 'location', 'EastOutside' );
%print(['figures/CorrMAP_' num2str(txt) '.ps'],'-depsc','-r500');

    %if  baz == 0 | baz <= 97 
    
    % for testing the quadrants %%%%%
% % % % % %     switch ( flag )
% % % % % % 			
% % % % % % % 			case 'residuals'
% % % % % %             case 'raw'
% % % % % %                 
% % % % % %           fid = fopen('ALL_raw_Q3RESIDS','a');
% % % % % %           fprintf(fid , ' %s  %f  %f  %f %f \n', sta, min(d), max(d), std(d), mean(d) );
% % % % % %          fclose(fid);
% % % % % %             case 'corrected'
% % % % % %                  fid = fopen('ALL_raw_Q3CORRECTED_RESIDS','a');
% % % % % %          fprintf(fid , ' %s  %f  %f  %f  %f \n', sta, min(d), max(d), std(d), mean(d) );
% % % % % %         fclose(fid);
% % % % % %         otherwise
% % % % % % 				
% % % % % % 				error( 'PTOMO:look_delays:bad_flag', ...
% % % % % % 					'flag must be ''corrected'', or ''residuals''' );
% % % % % %     end
    %end



end
% % % for iii=sta(1):4:sta(end)
% % % disp(['The plot is: ' num2str(sta(iii))])
% % %     %print(['figures/CorrMAP_' num2str(sta(iii)) '.pdf'],'-dpdf','-r600');
% % % end
