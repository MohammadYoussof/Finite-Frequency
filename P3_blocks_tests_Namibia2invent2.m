%%% Synthetic, STRUCTURE_DRIVEN blocks tests,  M.Y. April 2014 %%%


%load Namibia_P_Dataset_v03_W_normalVELOCITY
 load block_mod_invention1.mat
% mf = zeros( ny, nx, nz ); 
%  [modx,mody,modz] = meshgrid( par.modx, par.mody, par.modz );
% 
% 
% 
% % % % %%%%  A%   Northern    FAST    Anomaly #1 Namibia %%%EASTERN SIDE%%%%
%  xw=10; %ynw=40; 
%     
%     for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=300
%          
%          x1=xw;
%          x2=x1+200;%-(z*0.001);  %ANOMALY WIDTH
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=0.1*z+50+(xw-x)*(-1.015);
%             y2=ya+330; y1=y2-300;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=-amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
% 
%          end 
%     end
%     
% % %     
%     
%     
%     
% %%%   B%compliment to the NORTHERN FAST anomaly #1 Namibia   EASTERN ANOMALY
% % %   
%     xw=100; %ynw=300; 
%     
%     for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=300
%          
%          x1=xw-50 ;  %STARTING POSITION FOR X
%          x2=x1+450;%+z;
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=2.1*z+50+(xw-x)*(-1.0015);
%              y2=ya+0; y1=y2-750;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=-amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
% 
%          end 
%     end
% 
%     
%     
%     
%     
%   %%%  blockkkkkk    the TRANSITIONAL FAST Anomaly #
%      xw=-10; %ynw=100; 
% %     
%     for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=200
%          
%          x1=xw-350;
%          x2=x1+750;
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=-0.13*z-150+(xw-2*x)*0.70;
%             y2=ya+200; y1=y2-300;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=amp*0; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
%          end
%      end
% % % % 
% 
% 
% % %%%good blocckkkk
%      xw=-350; %ynw=100; 
%     
%     for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=200
%          
%          x1=xw-350;
%          x2=x1+750;
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=-0.13*z-150+(xw-2*x)*0.70;
%             y2=ya+200; y1=y2-600;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
%          end
%     end
%     
% %%%  blockkkkk  
% xw=-200; %yn
% 
%        for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=200
%          
%          x1=xw-250;
%          x2=x1+400;
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=0.13*z-150+(xw-2*x)*0.222;
%             y2=ya+500; y1=y2-400;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
%          end
%        end
% 
%      
%        
%        
%        
%        
%        xw=300; %yn
% 
%        for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=200
%          
%          x1=xw-250;
%          x2=x1+400;
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=0.13*z-150+(xw-2*x)*1.0111222;
%             y2=ya+480; y1=y2-550;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=-amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
%          end
%      end
%        
%        
%        
%        
%        
% 
% % % % %   %%%%%%%       WESTERN MIDDLE ANOMALY
%      xw=150; %ynw=-400;    
%      for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=300
%          
%          x1=xw-700; %x1 position
%          x2=x1+450; %-(z*0.0001); %x2 position
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%             ya=.269*z-350+(xw-x)*(0.2145);  %positioning y
%              
%             y2=ya+500; y1=y2-800;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
% 
%          end  
%      end
% % % %%%%%%%%ESTERN MIDDLE ANOMALY
%      
%       xw=-400; %ynw=-400;    
%      for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=300
%          
%          x1=xw+600; %x1 position
%          x2=x1+150; %-(z*0.0001); %x2 position
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%            ya=.1769*z-150+(xw-1*x)*(0.007272145);  %positioning y
%             %ya=0.13*z-150+(xw-2*x)*0.222; 
%             y2=ya+550; y1=y2-900;   
%             
%             
%          %   ya=-0.13*z-150+(xw-2*x)*0.70;
%          %   y2=ya+200; y1=y2-300;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=-amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
% 
%          end  
%      end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  xw=-200; %ynw=-400;    
%      for iz=1:nz
%          z=par.modz(iz);
%          
%          if z<=300
%          
%          x1=xw+600; %x1 position
%          x2=x1+150; %-(z*0.0001); %x2 position
%             
%          for ix=1:nx
%             if par.modx(ix)>x1 && par.modx(ix)<x2
%             x=par.modx(ix);
%            ya=.1769*z-150+(xw-1*x)*(0.007272145);  %positioning y
%             %ya=0.13*z-150+(xw-2*x)*0.222; 
%             y2=ya+550; y1=y2-900;   
%             
%             
%          %   ya=-0.13*z-150+(xw-2*x)*0.70;
%          %   y2=ya+200; y1=y2-300;   
%             
%             for iy=1:ny
%                 if par.mody(iy)>y1 && par.mody(iy)<y2
%                     mf(iy,ix,iz)=amp; 
%                     alerta=1;
%                 end    
%             end
% 
%             end
% 
%          end
% 
%          end  
%      end

% %  x=par.modx(ix);
% %             ya=-0.13*z-150+(xw-2*x)*0.70;
% %             y2=ya+200; y1=y2-300;  
 
    % if the output is m and GR and dR 
    % Then 
    
    %df = G*mf(:);
    new_df = GR*m(:);
    
    %noise = randn(size(df))*0.25;% Gaussian noise, standard deviation 
    %df = df + noise; % Add the noise
    DR = speye(nmodel);
    new_dR = [new_df; zeros(size(LR,1),1)];
    new_GR = [GR; DR*damp];
    new_dR = [new_dR; zeros(nmodel,1)];
    new_GR = [new_GR; LR*smooth];
    [new_m,flag,relres,iter,resvec] = lsqr( new_GR, new_dR, 1e-6, 100 );
    % ---- plot synthetic inversion ------------
    make_DS_map_synth_invent(m,new_m,par,stax,stay); 
    upsample_model_syn_m_invent(m,par,1);
    toc
    upsample_model_syn_new_m_invent(m,par,1);
    toc
    save block_mod_invention2