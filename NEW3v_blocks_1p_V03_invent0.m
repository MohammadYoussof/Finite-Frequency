% make vertical cross sections through S-tomography model(2012, JAN)
% make vertical cross sections through tomography model
clear all
clf
cmap=colormap(redblue);
cmap=colormap(flipud(colormap));
colormap(cmap);
load figures/plot_pts_syn_m_invent.mat

plot_pt_syn_m_invent.lat = plot_pts_syn_m_invent.lat';
plot_pt_syn_m_invent.lon = plot_pts_syn_m_invent.lon';
plot_pt_syn_m_invent.z = plot_pts_syn_m_invent.z';
plot_pt_syn_m_invent.val = plot_pts_syn_m_invent.val';

% pick endpoints of x section on map (give coordinates in lon,lat)

%A = [17.00 -34.25]; B=[29.50 -18.50]; %NORTH% A1
%A = [19.25 -34.25]; B=[31.50 -18.50];%MIDDLE% A2
%A = [21.50 -34.25]; B=[33.75 -18.50];%SOUTH% A3

%A = [8.00 -21.50]; B=[16.00 -18.00]; %General South_Africa% A


A = [7.00 -18.60]; B=[16.00 -18.50]; %  NAMIBIA PROFILE #1
%A = [17.00 -20.001]; B=[33.00 -20.00];
%A = [17.00 -27.501]; B=[33.00 -27.50];
for kk = 1:size(A,1)
    

    % subset input model around x-section
    plot_pt_syn_m_invent.dist1 = distance( plot_pt_syn_m_invent.lat, plot_pt_syn_m_invent.lon, ...
                    A(kk,2).*ones(size(plot_pt_syn_m_invent.lat)), A(kk,1).*ones(size(plot_pt_syn_m_invent.lat)) );
    plot_pt_syn_m_invent.dist2 = distance( plot_pt_syn_m_invent.lat, plot_pt_syn_m_invent.lon, ...
                    B(kk,2).*ones(size(plot_pt_syn_m_invent.lat)), B(kk,1).*ones(size(plot_pt_syn_m_invent.lat)) );
    xd = distance( A(kk,2),A(kk,1),B(kk,2),B(kk,1) );       % define an along Xsection distance
    xdr=(pi/180)*xd;
    plot_pt_syn_m_invent.dist12 = plot_pt_syn_m_invent.dist1 + plot_pt_syn_m_invent.dist2;
    sb_indx = find(plot_pt_syn_m_invent.dist12 < (xd+0.35) );
    xxx = plot_pt_syn_m_invent.lon(sb_indx);
    yyy = plot_pt_syn_m_invent.lat(sb_indx);
    zzz = -plot_pt_syn_m_invent.z(sb_indx);
    vp  = -plot_pt_syn_m_invent.val(sb_indx); % there's an issue here with the sign
    %---i changed this vp myslef...remember ya mohamed%%%%
    %%%vp  = plot_pt_syn_m_invent.val(sb_indx); % use this for preferred model, hte previous one for synth
    % depth range of section
    Z=[-1000 00];

    %build grid of points to interpolate model onto
    %h=100;%# horizontal increments (h+1=# of intervals)
    %v=50;%# vertical increments ""   ""
    %hi = 0.02;
    hi=0.12;
    v = round( (abs(Z(1))-abs(Z(2))) /12 );
    sect_length = distance(A(kk,2), A(kk,1), B(kk,2), B(kk,1));
    %sect_length = sqrt( (abs(A(kk,1))-abs(B(kk,1)))^2 + abs( A(kk,2)-B(kk,2))^2 );
    h = round(sect_length/hi);

    if A(kk,2)<B(kk,2) %for sections that are oriented SW-NE
        xx=A(kk,1):abs((A(kk,1)-B(kk,1)))/h:B(kk,1); xx=xx';%define lon points in map view(x)
        yy=A(kk,2):abs((A(kk,2)-B(kk,2)))/h:B(kk,2); yy=yy';%define lat points in map view(y)
        zz=Z(1):(Z(2)-Z(1))/v:Z(2);        %define points at depth (one column)
        % put 3 vectors above into 3 equal size matrices
        x=xx*ones(1,v+1);  %v+1 to match num of intervals
        y=yy*ones(1,v+1);
        z=ones(h+1,1)*zz;
    elseif A(kk,2)>B(kk,2) % NW-SE cross sections
        xx=A(kk,1):abs((A(kk,1)-B(kk,1)))/h:B(kk,1); xx=xx';
        yy=A(kk,2):-abs((A(kk,2)-B(kk,2)))/h:B(kk,2); yy=yy';
        zz=Z(1):(Z(2)-Z(1))/v:Z(2);        
        x=xx*ones(1,v+1);  
        y=yy*ones(1,v+1);
        z=ones(h+1,1)*zz;
    end

    %define an along Xsection distance
    hh = [ 0:sect_length/h:sect_length ]';
    %hh  = distance(y,x,A(kk,2).*ones(size(x)),A(kk,1).*ones(size(x)) );
    hhr = (pi/180).*hh;

    %interpolate velocity model onto points
     XVp = griddata(xxx,yyy,zzz,vp,x,y,z);
    % ----------------- plot -------------------------------------------------
    clf
   
    figure(1) 
    cmap=colormap(redblue);
    cmap=colormap(flipud(colormap));
    colormap(cmap);
    %box off;
    hhrp=((pi/2)-(pi/180)*0.5*xd)+hhr;
    r = [ 6371+(zz) ]';
    theta = hhrp(:,1)';
    X = -r*cos(theta);
    Y = r*sin(theta); 
    
    C = -XVp'.*100;
    contour(X,Y,C,500)
    
    
    cmap=colormap(redblue);
    cmap=colormap(flipud(colormap));
    set(gca,'clim', [-0.8 0.8])
    
    %pcolor(X,Y,C)
    %shading interp
    %box off;
    
    pcolor(X,Y,C)
    Csyn_bl=C;
    Xsyn_bl=X;
    Ysyn_bl=Y;
    
      colorbar
    shading flat
    axis equal
    axis tight
    hold on
    Y4 = (6371-0)*sin(hhrp(:,1));
    X4 = (6371-0)*cos(hhrp(:,1));
    plot(X4,Y4,'k--')
    Y6 = (6371-1000)*sin(hhrp(:,1));
    X6 = (6371-1000)*cos(hhrp(:,1));
    plot(X6,Y6,'k--')
    
%     Y8 = (6371-410)*sin(hhrp(:,1));
%     X8 = (6371-410)*cos(hhrp(:,1));
%     plot(X8,Y8,'r--')
%     Y10 = (6371-660)*sin(hhrp(:,1));
%     X10 = (6371-660)*cos(hhrp(:,1));
%     plot(X10,Y10,'r--')
    
% %     if zz < 6100 
% %     vvv = [0.9 0.9];
% %     contour(X,Y,C,vvv,'k--', 'LineWidth',1)
% %     end
% %     
% %     if zz < 6100 
% %     vvv = [0.75 0.75];
% %     contour(X,Y,C,vvv,'r--', 'LineWidth',1)
% %     end
% %     
    %I added this bit, don't know if it works%%%good
     yt=[5371:100:6371];
     
    set(gca,'YTick',yt,'YTickLabel',(6371-yt))
    xt=[-500:100:500];
    set(gca,'XTick',xt,'XTickLabel',(xt))
    axis([-500 500 5371 6371])
    
     colormap(redblue)
    colormap(flipud(colormap))
    set(gca,'clim', [-0.8 0.8])
 
    cs = max(abs([min(C(:)),max(C(:))]));
     caxis([-cs cs])
set(gcf,'units',get(gcf,'PaperUnits'),'Position',get(gcf,'PaperPosition'))
    set(gca,'clim', [-0.8 0.8])
  
     print(['figures/P_blocks_m_model_syntha.pdf'],'-dpdf','-r950');
     
    
end
%%%%save blocks_synth1