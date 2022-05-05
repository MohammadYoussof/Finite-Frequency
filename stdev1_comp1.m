load observ.mat
load synth.mat


load figures/plot_pts.mat

plot_pt.lat = plot_pts.lat' ;
plot_pt.lon = plot_pts.lon' ;
plot_pt.z = plot_pts.z' ;
plot_pt.val = plot_pts.val';
figure(1)

DX1=[13.00 -16.50]; DX2 =[15.75 -22.0]; %%%%%%FF'
%DX1=[11.00 -16.50]; DX2 =[13.75 -22.0]; %%%%EE'%%%%
%DX1=[8.00 -16.50]; DX2 =[11.75 -22.0]; %%%%DD''%%%%

%%good alternative%%DX1=[8.50 -16.00]; DX2 =[11.75 -22.50]; %%%%dd''%%%%


%DX1=[7.50 -17.00]; DX2 =[10.00 -23.00];  %NS


%%%%DX1 =[13.00 -16.00]; DX2=[15.75 -23.00];   %// coast


 %simialr to A:: :: %%DX1=[8.00 -20.00]; DX2 =[15.00 -17.50]; % XV01A= CROSSING CONT-OCEANIC TRANSITION!

 %need reversal %DX1 =[12.50 -16.00]; DX2=[15.75 -23.00];   %// coast= NE-SW

 %&&&&&&&&&&&&  SELECTED  PROFILES#####

%%%%%%SW-NE profiles through the oceanic - continent // to the WR %%%%%
% % DX1 =[7.00 -23.00]; DX2=[16.00 -19.00];   %// C %%%%
% % 
% % DX1 =[7.00 -22.00]; DX2=[16.00 -18.00];  %%%  B %%%
% % 
%DX1=[7.00 -20.00]; DX2 =[15.00 -17.00]; % %%A   %%%%XV01= CROSSING CONT-OCEANIC TRANSITION!


for kk = 1:size(DX1,1)
    

    % subset input model around x-section
    plot_pt.dist1 = distance( plot_pt.lat, plot_pt.lon, ...
                    DX1(kk,2).*ones(size(plot_pt.lat)), DX1(kk,1).*ones(size(plot_pt.lat)) );
    plot_pt.dist2 = distance( plot_pt.lat, plot_pt.lon, ...
                    DX2(kk,2).*ones(size(plot_pt.lat)), DX2(kk,1).*ones(size(plot_pt.lat)) );
    xd = distance( DX1(kk,2),DX1(kk,1),DX2(kk,2),DX2(kk,1) );       % define an along Xsection distance
    xdr=(pi/180)*xd;
    plot_pt.dist12 = plot_pt.dist1 + plot_pt.dist2;
    sb_indx = find(plot_pt.dist12 < (xd+0.35) );
    xxx = plot_pt.lon(sb_indx);
    yyy = plot_pt.lat(sb_indx);
    zzz = -plot_pt.z(sb_indx);
    %vp  = -plot_pt.val(sb_indx); % there's an issue here with the sign
    %---i changed this vp myslef...remember ya mohamed%%%%
    vp  = plot_pt.val(sb_indx); % use this for preferred model, the previous one for synth
    % depth range of section
    Z=[-350 00];

    %build grid of points to interpolate model onto
    %hi = 0.02;
    hi=0.12;
    v = round( (abs(Z(1))-abs(Z(2))) /12 );   %# vertical increments ""   ""
    sect_length = distance(DX1(kk,2), DX1(kk,1), DX2(kk,2), DX2(kk,1));
    %sect_length = sqrt( (abs(DX1(kk,1))-abs(DX2(kk,1)))^2 + abs( DX1(kk,2)-DX2(kk,2))^2 );
    h = round(sect_length/hi);   %horizontal increments (h+1=# of intervals)

    if DX1(kk,2)<DX2(kk,2) %for sections that are oriented SW-NE
        xx=DX1(kk,1):abs((DX1(kk,1)-DX2(kk,1)))/h:DX2(kk,1); xx=xx';%define lon points in map view(x)
        yy=DX1(kk,2):abs((DX1(kk,2)-DX2(kk,2)))/h:DX2(kk,2); yy=yy';%define lat points in map view(y)
        zz=Z(1):(Z(2)-Z(1))/v:Z(2);        %define points at depth (one column)
        % put 3 vectors above into 3 equal size matrices
        x=xx*ones(1,v+1);  %v+1 to match num of intervals
        y=yy*ones(1,v+1);
        z=ones(h+1,1)*zz;
    elseif DX1(kk,2)>DX2(kk,2) % NW-SE cross sections
        xx=DX1(kk,1):abs((DX1(kk,1)-DX2(kk,1)))/h:DX2(kk,1); xx=xx';
        yy=DX1(kk,2):-abs((DX1(kk,2)-DX2(kk,2)))/h:DX2(kk,2); yy=yy';
        zz=Z(1):(Z(2)-Z(1))/v:Z(2);        
        x=xx*ones(1,v+1);  
        y=yy*ones(1,v+1);
        z=ones(h+1,1)*zz;
    end

    %define an along Xsection distance
    hh = [ 0:sect_length/h:sect_length ]';
    %hh  = distance(y,x,DX1(kk,2).*ones(size(x)),DX1(kk,1).*ones(size(x)) );
    hhr = (pi/180).*hh;
    
    %interpolate velocity model onto points
     XVp = griddata(xxx,yyy,zzz,vp,x,y,z);
    % ----------------- plot -------------------------------------------------
  
     %figure(2) 
   
    subplot(2,1,1)

    [sname,slat,slon,sz]=textread('sta1.file','%s %f %f %f');
    plot(slon,slat,'k.')
    hold on
    


    plot(([DX1(kk,1) DX2(kk,1)]),([DX1(kk,2) DX2(kk,2)]),'r','Linewidth',3)
    axis equal, axis manual
     axis( [ 7 18 -23 -16 ] )
   
    
    coast = load('coastlines');
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
    title( sprintf('Stations map and selected profiles positions'))
  
    hold on
       
    hb=subplot(2,1,2);
    
    
    %%%---
    %figure(1)
    
    hhrp=((pi/2)-(pi/180)*0.5*xd)+hhr;
    
    r = [ 6371+zz ]';
    theta= hhrp(:,1)';
    X = -r*cos(theta);
    Y = r*sin(theta); 
    C = -XVp';
    
    
    
    Xn=X;
    Yn=Y;
    Cn=Cd-Co;
    
    %Csd=stdfilt(Cn);
    
    %Cn=Cn*0;
    
    %Cn=Cn+Csd;
    
    contourf(Xn,Yn,Cn,500)
  
   % set(gca,'clim', [min(plot_pts.val) max(plot_pts.val)])
    set(gca,'clim', [-0.8 0.8])
    pcolor(Xn,Yn,Cn)
     colorbar
    %shading interp
    shading flat
    axis equal
    axis tight
    hold on
    Y0 = (6371-0)*sin(hhrp(:,1));
    X0 = (6371-0)*cos(hhrp(:,1));
    plot(X0,Y0,'-.m','LineWidth',0.5)
    Y1 = (6371-100)*sin(hhrp(:,1));
    X1 = (6371-100)*cos(hhrp(:,1));
    plot(X1,Y1,'-.m','LineWidth',0.5)
    Y2 = (6371-200)*sin(hhrp(:,1));
    X2 = (6371-200)*cos(hhrp(:,1));
    plot(X2,Y2,'-.m','LineWidth',0.5)
    Y3 = (6371-300)*sin(hhrp(:,1));
    X3 = (6371-300)*cos(hhrp(:,1));
    plot(X3,Y3,'-.m','LineWidth',0.5)
    Y4 = (6371-400)*sin(hhrp(:,1));
    X4 = (6371-400)*cos(hhrp(:,1));
    plot(X4,Y4,'-.m','LineWidth',0.5)
    Y5 = (6371-500)*sin(hhrp(:,1));
    X5 = (6371-500)*cos(hhrp(:,1));
    plot(X5,Y5,'-.m','LineWidth',0.5)

    
    
    %I added this bit, don't know if it works%%%good
    yt=[6021:100:6371];
    set(gca,'YTick',yt,'YTickLabel',(6371-yt))
    %set(gca, 'YTick', []);
     xt=[-390:100:390];
    set(gca,'XTick',xt,'XTickLabel',(xt))
    
    
%     ss=[X(700) Y(700)];
%     ww=[X(700) Y(1)];
%     plot(ss,ww,'k--','linestyle','-')
    %axis([-340 340 5630 6371])  %sample
    axis([-390 390 6021 6371])
    %cbfit(5,0)
    colormap(redblue)
    colormap(flipud(colormap))
    set(gca,'clim', [-0.8 0.8])
   cs = max(abs([-0.8, 0.8]));
    %cs = max(abs([min(C(:)),max(C(:))]));
     caxis([-cs cs])
    set(gcf,'units',get(gcf,'PaperUnits'),'Position',get(gcf,'PaperPosition'))
     %xlabel( 'Profiles Length in km, 0 is (lat -26.88 & lon 25.325)' );
     
     title({'Depth slices through the P-wave velocity model'; 'Ocean Profile'})
     set(hb, 'xscale', 'linear');
    print('figures/final_VXsections/v01A_comp1SD1okokoko.pdf','-dpdf','-r950');
   
end
%save stdcomp2
