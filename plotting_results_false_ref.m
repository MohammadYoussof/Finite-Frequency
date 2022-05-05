% plot the true model results ##########


clear 
close all
clc

disp('AsSalam Alykom Dr. Mohammad bin Youssof :) :) :) ')  
disp(' This is for PLOTTING Tomography Inversion RESULTS')
disp(' --------------------------------------------')
 

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

%==========***** subfunction options/inversion parameters*****==========%%%
plot_delays    = true; %Plotting original delaytimes%%
site_corr      = true; % calculate crust corrections%%
uniform_grid   = false; %Gridding style%%
plot_nodes     = true; %Plotting nodes spacings%%
plot_hitq      = false;%Plotting hitting quality %%
build_G        = true; %Building G Matrix %%
synth_test     = false; %Implementing Resolution Tests %
build_LR       = true; %Building Regularization Matrix %%
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





load figures_false_ref/Namib_False_Ref.mat

figure(1)

if plot_tomo_maps == true
   
    make_DS_map_false(vm,par,stax,stay);
end

%tic

%toc


figure(2)
 
%%%


%         %plot station residuals%%
% 

  
 s_residuals_false=stn.resid;
 plot_s_resids_false(par, data, s_residuals_false);  hold on; set(gca,'clim', [-0.5 0.5])
   cs = max(abs([-0.8, 0.8]));
  print(['figures_false_ref/stns_residuals_false.pdf'],'-dpdf','-r500');

 
% load recObj123
% play(recObj)


