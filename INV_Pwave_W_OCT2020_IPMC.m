%%%%%%%%%%%%%%%%%%   Main Code of the P-WAVE Tomography%%%%%%%%%%%%%%%%%%%
%-----------------my work in this code started on June 2010----------%%%%
%%%%%%%%%% Max Bezada ----  Thank you for teaching me Tomography ----%%%%
%-----------------last update :  August, 20th , 2013-----------------%%%%
%-----------------last update :  Spet, 21th , 2015-----------------%%%%%
%-----------------last update :  Oct, 20th , 2020-----------------%%%%%

clear 
close all
clc

disp('AsSalam Alykom Dr. Mohammad bin Youssof :) :) :) ')  
disp(' This is the P-WAVE Tomography Inversion Code')
disp(' --------------------------------------------')
disp(' WALPASS Land and OBS Dataset')

%----- Define geographic model space and tomographic method parameters%%%

origin = [ -19.0 12 ];          % Array center [ lat lon ]
zmax   = 1005;                        % Max depth in model domain was 1200 %
mstruct = defaultm('mercator');      % Projection%
mstruct.origin = [origin 0];
mstruct = defaultm( mstruct );
par = struct( ...
'verbose', 2, ...
'map_proj', mstruct, ...
'origin', origin );


%--------------------color scheme----------------------------------%%%
cmap = colormap(redblue);
cmap = flipud(cmap);
colormap(cmap);
%-------------------------------------------------------------------%%%
load Namibia_P_Dataset_ref1
%==========***** subfunction options/inversion parameters*****==========%%%
plot_delays    = true; %Plotting original delaytimes%%
site_corr      = true; % calculate crust corrections%%
uniform_grid   = false; %Gridding style%%
plot_nodes     = true; %Plotting nodes spacings%%
plot_hitq      = false;%Plotting hitting quality %%
build_G        = false; %Building G Matrix %%
synth_test     = false; %Implementing Resolution Tests %
build_LR       = false; %Building Regularization Matrix %%
output_model   = true; %Calculating the output model %%
damp           = 2; %42; %%%%%16;   %Damping Factor %% %was 5
smooth         = 10; %49; %%%%%25;   %Normal Damping Factor %% %was 8,,,,,,15
par.damp_stn   = 10.6;  %If you want stn statics off, set this >>1 %%
par.damp_evt   = 0.0;  %I think if you want event_statics off, set this >>1 

par.use_crust_corr = true; % apply crust corrections in inversion

plot_tomo_maps = true;     % Plotting final models %%
plot_tomo_syn_maps = true; % Plotting resolution results %%
plot_stn_static= true;     % Plotiing station static terms %%
plot_resids    = true;     % Plotting final residuals %%

%%%%%==========***** "Building G Matrix" *****============%%%%%
%if build_G == true
    %=================== Input Rays ======================================%
    %fid = fopen('walpass2.file', 'r');
    fid = fopen('n3.dat', 'r');
   %%%%%% fid = fopen('offshore_sea.file', 'r');
    A = [textscan(fid, '%f %f %f %s %d %f %f %f %f %d'),1]; fclose(fid);
    ray.pd    = A{1};
    ray.p     = ray.pd .*(360/(2*pi*6371));
    ray.baz   = A{2};
    %ray.d     = A{3};
    dc1      =dc;
    ray.d    = dc1;
    save dc1 dc1
    
    %NICE OMO
    for i = 1:length(ray.baz)
        if ray.baz(i)> 180 &&  ray.baz(i) <= 270
            ray.d(i)=ray.d(i)*0.5;
        end
    end  
    ray.sta   = A{4};
    ray.orid  = A{5};
    ray.cf    = A{6};
    ray.chan  = A{10};
    ray.lat   = A{7};
    ray.lon   = A{8};
    sta_ray   = ray.sta;
    ray.nrays = length(ray.p);
    
% %     ray.pd(ray.baz<180)= [];
% %     ray.p(ray.baz<180)= [];
% %     ray.d(ray.baz<180)= [];
% %     ray.sta(ray.baz<180)= [];
% %     ray.orid(ray.baz<180)= [];
% %     ray.cf(ray.baz<180)= [];
% %     ray.chan(ray.baz<180)= [];
% %     ray.lat(ray.baz<180)= [];
% %     ray.lon(ray.baz<180)= [];
% %     sta_ray   = ray.sta;
% %     ray.nrays = length(ray.p);
% %     ray.baz(ray.baz<180)= [];
% %     
% % %%%ray.pd(ray.baz<180)= [];
% % %%%ray.p(ray.baz<180)= [];
% % %%%ray.d(ray.baz<180)= [];
% % %%%ray.sta(ray.baz<180)= [];
% % %%%ray.orid(ray.baz<180)= [];
% % %%%ray.cf(ray.baz<180)= [];
% % %%%ray.chan(ray.baz<180)= [];
% % %%%ray.lat(ray.baz<180)= [];
% % %%%ray.lon(ray.baz<180)= [];
% % %%%%sta_ray   = ray.sta;
% % %%%ray.nrays = length(ray.p);
% % %%%ray.baz(ray.baz<180)= [];

    
    %=================== Input Stations ==================================% 
    sta_names = unique( ray.sta );
% % % % %    %fid = fopen('new_crust_model_modOCEAN.txt','r'); % THICKER crust%TEST%%%%
% % % % %    %fid = fopen('new_crust_model_modOCEAN.txt','r');  %THINNER crust % TEST-2%%%
% % % % %    
% % % % %    %fid = fopen('thicknew.dat','r');
% % % % %     
    fid = fopen('new_crust_model.txt','r'); %correct%%
    %fid = fopen('new_crust_model555.txt','r'); %correct%% T=JUST FOR CHECKING ..TEMPPPPPPP
    %fid = fopen('new_crust_model_THINOCEAN.txt','r');   %correcr+5v%
    %fid = fopen('new_crust_model_THICKOCEAN.txt','r');   %correcr-5v%
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
    
      %------ Station Corrections-------------------%
%  if site_corr==true
%         
%     data = crust_corr(data);
%  end 

    ray = data.ray;
    stn = data.stn;
% else
%     load data data;
%     ray = data.ray;
%     stn = data.stn;
%     node = data.node;
% end
            %--------%Plotting original delaytimes, after crustal corrections----------------%
if plot_delays == true
  stds = look_delays( par, data, 'corrected' );
end
save stds stds

%%%%%==***** "Parameterization of Nodes Spacings"*****===========%%%%
%--------------define model space ------------------------%
    [stax,stay] = project_xy( par, stn.lat, stn.lon );
    data.stn.stax = stax;
    data.stn.stay = stay;
    
  if uniform_grid==true
        dh = 30;
        margin1 = 100;
        xmin = floor( min(stax)/dh )*dh-margin1;
        xmax = ceil( max(stax)/dh )*dh+margin1;
        ymin = floor( min(stay)/dh )*dh-margin1;
        ymax = ceil( max(stay)/dh )*dh+margin1;
        modx = xmin:dh:xmax;
        mody = ymax:-dh:ymin;
        modz = [ 1 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000]; % ok for 700 km domain%
        if zmax>max(modz)
            zmax = max(modz);
        end
        modz = modz(modz<=zmax);
        modz = [ modz max(modz)+10 ];
 %m
        nx = length( modx );
        ny = length( mody );
        nz = length( modz );
        nmodel = nx*ny*nz;

        par.modx = modx;
        par.mody = mody;
        par.modz = modz;
        par.nx = nx;    par.ny = ny;    par.nz = nz;
        par.nmodel = nmodel;
%m
  elseif uniform_grid==false
        dh1 = 30;  %35
        dh2 = 35;  %40
        dh3 = 40;  %45
        dh4 = 45;  %50
        margin2 = 50; %it was 150
        modx1 = round(min(stax))+margin2:dh1:round(max(stax))-margin2;
        mody1 = round(max(stay))-margin2:-dh1:round(min(stay))+margin2;
        modx2 = [ min(modx1)-2*dh2 min(modx1)-dh2 modx1 max(modx1)+dh2 max(modx1)+2*dh2 ];
        mody2 = [ max(mody1)+2*dh2 max(mody1)+dh2 mody1 min(mody1)-dh2 min(mody1)-2*dh2 ];
        modx3 = [ min(modx2)-dh3 modx2 max(modx2)+dh3 ];
        mody3 = [ max(mody2)+dh3 mody2 min(mody2)-dh3 ];
        xmin = min(modx3)-3*dh4;
        xmax = max(modx3)+3*dh4;
        ymin = min(mody3)-3*dh4;
        ymax = max(mody3)+3*dh4;
        modx = [ xmin:dh4:min(modx3) modx3(2:end-1) max(modx3):dh4:xmax];
        mody = [ ymax:-dh4:max(mody3) mody3(2:end-1) min(mody3):-dh4:ymin ];
  	    %modz = [ 1 15 36 50 85 100 125 150 175 200 225 250 275 300 325 350 375 400 450 500 550 600 700 800 900 1000]; % ok for 700 km domain%
        modz = [ 1 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000]; % ok for 700 km domain%
       
        if zmax>max(modz)
            zmax = max(modz);
        end
        modz = modz(modz<=zmax);
        modz = [ modz max(modz)+10 ];

        nx = length( modx );
        ny = length( mody );
        nz = length( modz );
        nmodel = nx*ny*nz;

        par.modx = modx;
        par.mody = mody;
        par.modz = modz;
        par.nx = nx;    par.ny = ny;    par.nz = nz;
        par.nmodel = nmodel;
  end
            %--------%Plotting nodes spacings-----------%
    if plot_nodes==true
        clf
        map_nodes(par, stax, stay);

    end    
%%%%%==========***** "Completing the G Matrix" *****=================%%%%%
if build_G == false
    load G G;   
elseif build_G == true
   [G,par,data] = make_G( par, data ); 
   data = node_hit_q(data);
   node = data.node;
   save data data;
   save G G;
end            
            % ----------%plot node hit quality -----------%
% if plot_hitq == true
% make_DS_map_hitq(node.hit_q,par,stax,stay);
% end

if plot_hitq == true
make_DS_map_hitq_real_rays(node.hit_q,par,stax,stay);
end
%%%%%==========***** "Regularization Martix"*****=================%%%%%
if build_LR==true
    LR = reg_mat(par);
    save LR LR;
elseif build_LR==false
    load LR LR;
end

%%%%%==========***** "Inverting synthetics and/or data" *****=============%%%%%
            %----------- resolution test --------%
           % ---- making checkerboard pattern ----- %
if synth_test == true
    mf = zeros( ny, nx, nz ); 
    [modx,mody,modz] = meshgrid( par.modx, par.mody, par.modz );

    amp = 0.008; 
    width = 4;
    rev_dep = 5;
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
    save mf_syn mf;

                % ---- invert synthetic data ------%
    df = G*mf(:);
    
    noise = randn(size(df))*0.15; %gaussian noise , sd randn gaussain 
    
    df = df +noise;
    DR = speye(nmodel);
    dR = [df; zeros(size(LR,1),1)];
    GR = [G; DR*damp];
    dR = [dR; zeros(nmodel,1)];
    GR = [GR; LR*smooth];
    [m,flag,relres,iter,resvec] = lsqr( GR, dR, 1e-6, 100 );

    % ---- Plotting synthetic inversion ------%%
    make_DS_map_synth(mf,m,par,stax,stay);  
%elseif synth_test == false
 % load mf mf
%    load m m 
   
end
tic
%output SYNTH high-res pattern%
if synth_test == true
upsample_model_syn_m(m,par,1);
toc
tic

%output SYNTH high-res model%
upsample_model_syn_mf(mf,par,1);
toc
end

            % ----  Inverting Real Data ------%
num_iter = 150;
mean_corr_tt=mean(data.ray.d);
if par.use_crust_corr == true
    d=(data.ray.d - data.ray.corr1 - mean_corr_tt);
else
    d=(data.ray.d - mean_corr_tt);
end
norid = length(unique(data.ray.orid));
nstas = length(unique(data.ray.sta_num));   end_s_stat = nmodel+nstas;
nmat  = nmodel+nstas+norid;

%clear GR dR;  NOT SURE ,, LET'S TRY IT
[GR,dR] = regularized_matrix( G, data, par, damp, ...
  smooth, LR );
[m,flag,relres,iter,resvec] = lsqr( GR, dR, 1e-6, num_iter );
fprintf( 'the LSQR terminated iteration at %d of %d.\n', iter, num_iter );
s_static  = m(nmodel+1:end_s_stat);
e_static  = m(end_s_stat+1:end);
vm = m(1:nmodel);
df = G*vm(:);


%noise = randn(size(df))*2.5; %gaussian noise , sd randn gaussain  
%df = df +noise;


dc = calc_dc(s_static,e_static,df(:),ray,stn); 
vr    = variance_reduction( d, dc(:) );
fprintf( '\tVariance Reduction = %.3f%%\n', vr );

norm_m = norm( vm(:) );
resid = norm( dc(:) - d );
data.ray.resids = dc(:) - d;
r_resids = data.ray.resids; %did it M%%
RMS_DATA=rms(data.ray.d);
RMS_RESIDUALS=rms(r_resids);
% save s_static s_static;
% save e_static e_static;
% save vm vm;
% save vr vr;
% save df_real df ;
% save dc_real dc ;
% save norm_m norm_m;
% save resid_real resid ;
% save r_resids r_resids;
% save RMS_DATA RMS_DATA;
% save RMS_RESIDUALS RMS_RESIDUALS;
%%%%%==========***** "Output Model" *****=============%%%%%%
    %------- output lat, lon, depth, Dv text file ------%%%
if output_model==true
    fid = fopen('latlonr_model','w');
    n=1;
    for kk = 1:nz
        for jj = 1:nx
            for ii = 1:ny
               n_indx = fast_sub2ind3d( [ny nx nz], ii, jj, kk );
               val    = -100*vm(n_indx);
               xx     = modx(jj);
               yy     = mody(ii);
               zz     = modz(kk);
               r      = zz;
               [lt,ln]= project_xy(par, xx, yy, 'inverse');
               fprintf(fid,' %f  %f  %f  %f \n', lt,ln,r,val);
               n=n+1;
            end
        end
    end
    fclose(fid);
end

%%%%%==========***** " Display Results " ******========%%%%%
            %--- plot Tomographic map slices
if plot_tomo_maps == true
    figure(1);
    make_DS_map(vm,par,stax,stay);
end

tic
            %---output High-Resolution Model
upsample_model(vm,par,1);
toc


figure(2)
            %-- plot STATION residuals at each station separately%% 
if plot_resids == true
  stds = look_delays( par, data, 'residuals' );
end
%%%


        %plot station residuals%%

stn.resid = zeros(length(stn.num),1);
for ii = 1:length(stn.num)
    indx2 = find(ray.sta_num == stn.num(ii));
    stn.resid(ii) = sum( (dc(indx2)-d(indx2)) )/length(indx2); %#
end   
clf
figure(3)
 s_residuals=stn.resid;
 plot_s_resids(par, data, s_residuals);
% % %  %M% print(['figures/Station_residuals.pdf'],'-dpdf','-r400');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %plot station residuals   TRIAL%%

stn.resid = zeros(length(stn.num),1);
for ii = 1:length(stn.num)
    indx2 = find(ray.sta_num == stn.num(ii));
    stn.resid(ii) = sum( d(indx2) )/length(indx2); %# TRIALS
end   

figure(4)
 s_residuals=stn.resid;
 plot_s_resids(par, data, s_residuals);
 %M% print(['figures/Station_residuals_TRIAL1.pdf'],'-dpdf','-r400');
 
 
 %%%%%%%%%%%%%%%%%%%%%
 
 
            % plot station static terms%
if par.use_crust_corr == true
 plot_s_static(par, data, s_static);
end
figure(5);
%M% print(['figures/Station_Static.pdf'],'-dpdf','-r400');
            % plot STATION static histogram 

figure(6);
hist(s_static(:))
title('Station Static Terms')
xlabel('seconds')
%M% print(['figures/Station_Static_Histogram.pdf'],'-dpdf','-r400');

            % plot EVENT static histogram 
figure(7);
hist(e_static(:))
title('Event Static Terms')
xlabel('seconds')
%M% print(['figures/Svent_Static_Histogram.pdf'],'-dpdf','-r400');
save Namib_False_Ref
% load recObj123
% play(recObj)


