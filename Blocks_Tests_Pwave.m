load Pmodel_13thOCT

mf = zeros( ny, nx, nz ); 
 [modx,mody,modz] = meshgrid( par.modx, par.mody, par.modz );
%%% Synthetic M.Y. %%%
%synthVDH

% % % %%%%  A%   Anomaly #1 kaapvaal %%%WESTERN SIDE%%%%
 xw=-300; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<=300
         
         x1=xw;
         x2=x1+350;%-(z*0.001);  %ANOMALY WIDTH
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=0.1*z+50+(xw-x)*(-.015);
            y2=ya+300; y1=y2-500;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
    
   % hold on
%     %%%   B%compliment to anomaly #1 kaapvaal  BIG EASTERN ANOMALY
% % %     
    xw=-300; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<=300
         
         x1=xw+250 ;  %STARTING POSITION FOR X
         x2=x1+550;%+z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=0.1*z+50+(xw-x)*(-.0000015);
             y2=ya+450; y1=y2-750;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
%%%%%%  C%
    xw=150; ynw=200; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<=300
         
         x1=xw-200;
         x2=x1+750;%-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.3*z+500+(xw-x)*(-.145);
            %ya=-100+(xw-x)*-3;
            y2=ya+300; y1=y2-500;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
      %%%    D% upper right positive (fast) anomaly%%%
     xw=-50; ynw=380; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<180
         
         x1=xw;
         x2=x1+750-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.3*z+100+(xw-x)*(-.145);
            y2=ya+400; y1=y2-100;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
%     
%     
%     
%     %%%%% E %%%%ZIMBABWE ANOMALYYYYY%%%5

    
         xw=-100; ynw=200; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<=250
         
         x1=xw+300;
         x2=x1+550-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.3*z+650+(xw-x)*.15;
            y2=ya+350; y1=y2-150;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
% %     hold on
% %     
%   
% %     
% 
  %%%   F%  Southern kaapvaal anomaly%%%
     xw=-100; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z<=300
         
         x1=xw-350;
         x2=x1+550+z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.13*z-200+(xw-2*x)*.15;
            y2=ya+100; y1=y2-350;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end

% % %      
% % %      hold on
% % %   %%%%%%%% G% Namaqua Natal and CapeFold Belts
     xw=150; ynw=-400;    
     for iz=1:nz
         z=par.modz(iz);
         
         if z<=300
         
         x1=xw-700; %x1 position
         x2=x1+550; %-(z*0.0001); %x2 position
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=.269*z-450+(xw-x)*(0.2145);  %positioning y
            %ya=-100+(xw-x)*-3;
            y2=ya-00; y1=y2-600;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
     end
% %     

%    %%%    H%  Simple anomaly to cover the complicated signal at the right conrner
% 
    xw=150; ynw=600;    
     for iz=1:nz
         z=par.modz(iz);
         
         if z>180 && z<=300
         
         x1=xw+50; %x1 position
         x2=x1+550;%-z; %x2 position
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=.0003*z+650+(xw-x)*(.000145);  %positioning y
            %and 0.0003 is the movement rate in y direction 
            y2=ya+200; y1=y2-500;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=amp; 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
     end
     
   %%%%%%% #####  white anomalies except Kaapvaal block %%%%%%%
   % #####  white anomalies except Kaapvaal block %%%
   
   
   
% % % % % % % % % % % %    
% % % % % % % % % % % %    % % % %%%%A!   Anomaly #1 kaapvaal %%%WESTERN SIDE%%%%
% % % % % % % % % % % %  xw=-300; ynw=400; 
% % % % % % % % % % % %     
% % % % % % % % % % % %     for iz=1:nz
% % % % % % % % % % % %          z=par.modz(iz);
% % % % % % % % % % % %          
% % % % % % % % % % % %          if z>300
% % % % % % % % % % % %          
% % % % % % % % % % % %          x1=xw;
% % % % % % % % % % % %          x2=x1+350;%-(z*0.001);  %ANOMALY WIDTH
% % % % % % % % % % % %             
% % % % % % % % % % % %          for ix=1:nx
% % % % % % % % % % % %             if par.modx(ix)>x1 && par.modx(ix)<x2
% % % % % % % % % % % %             x=par.modx(ix);
% % % % % % % % % % % %             ya=0.1*z+50+(xw-x)*(-.015);
% % % % % % % % % % % %             y2=ya+300; y1=y2-500;   
% % % % % % % % % % % %             
% % % % % % % % % % % %             for iy=1:ny
% % % % % % % % % % % %                 if par.mody(iy)>y1 && par.mody(iy)<y2
% % % % % % % % % % % %                     mf(iy,ix,iz)=(amp); 
% % % % % % % % % % % %                     alerta=1;
% % % % % % % % % % % %                 end    
% % % % % % % % % % % %             end
% % % % % % % % % % % % 
% % % % % % % % % % % %             end
% % % % % % % % % % % % 
% % % % % % % % % % % %          end
% % % % % % % % % % % % 
% % % % % % % % % % % %          end %(end if z <250)
% % % % % % % % % % % %     end
% % % % % % % % % % % %     
% % % % % % % % % % % %    % hold on
% % % % % % % % % % % % %     %%%B! compliment to anomaly #1 kaapvaal  BIG EASTERN ANOMALY
% % % % % % % % % % % % % % %     
% % % % % % % % % % % %     xw=-300; ynw=400; 
% % % % % % % % % % % %     
% % % % % % % % % % % %     for iz=1:nz
% % % % % % % % % % % %          z=par.modz(iz);
% % % % % % % % % % % %          
% % % % % % % % % % % %          if z>300 
% % % % % % % % % % % %          
% % % % % % % % % % % %          x1=xw+250 ;  %STARTING POSITION FOR X
% % % % % % % % % % % %          x2=x1+550;%+z;
% % % % % % % % % % % %             
% % % % % % % % % % % %          for ix=1:nx
% % % % % % % % % % % %             if par.modx(ix)>x1 && par.modx(ix)<x2
% % % % % % % % % % % %             x=par.modx(ix);
% % % % % % % % % % % %             ya=0.1*z+50+(xw-x)*(-.0015);
% % % % % % % % % % % %              y2=ya+450; y1=y2-850;   
% % % % % % % % % % % %             
% % % % % % % % % % % %             for iy=1:ny
% % % % % % % % % % % %                 if par.mody(iy)>y1 && par.mody(iy)<y2
% % % % % % % % % % % %                     mf(iy,ix,iz)=(amp); 
% % % % % % % % % % % %                     alerta=1;
% % % % % % % % % % % %                 end    
% % % % % % % % % % % %             end
% % % % % % % % % % % % 
% % % % % % % % % % % %             end
% % % % % % % % % % % % 
% % % % % % % % % % % %          end
% % % % % % % % % % % % 
% % % % % % % % % % % %          end %(end if z <250)
% % % % % % % % % % % %     end


% % % %%%%  A%   Anomaly #1 kaapvaal %%%WESTERN SIDE%%%%
 xw=-300; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw;
         x2=x1+350;%-(z*0.001);  %ANOMALY WIDTH
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=0.1*z+50+(xw-x)*(-.015);
            y2=ya+300; y1=y2-600;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp*(-1); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
    
   % hold on
%     %%%   B%compliment to anomaly #1 kaapvaal  BIG EASTERN ANOMALY
% % %     
    xw=-300; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw+250 ;  %STARTING POSITION FOR X
         x2=x1+550;%+z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=0.1*z+50+(xw-x)*(-.0000015);
             y2=ya+450; y1=y2-750;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp*(-1); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end




  %%%%%%C!   %%%%
   xw=150; ynw=200; 
   
    for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw-200;
         x2=x1+750;%-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.13*z+500+(xw-x)*(-.145);  % changed
            %ya=-100+(xw-x)*-3;
            y2=ya+300; y1=y2-500;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=(amp*0); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
% % %    D!   %%%upper right positive anomaly%%%
     xw=-50; ynw=380; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw;
         x2=x1+750-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.3*z+100+(xw-x)*(-.145);
            %ya=-100+(xw-x)*-3;
            y2=ya+400; y1=y2-100;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=(-amp*0); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
    
%     
%     
% % % %     %%%E!  %% %%%%ZIMBABWE ANOMALYYYYY%%%5
% % % 
    
         xw=-100; ynw=200; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw+300;
         x2=x1+550-z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.3*z+650+(xw-x)*.15;
            y2=ya+350; y1=y2-150;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=(amp*0); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end
% %     hold on
% %     
%   
% %     
% 
% % % % % % % % % % %   %%%F!  %%%southern kaapvaal anomaly%%%
% % % % % % % % % % %      xw=-100; ynw=400; 
% % % % % % % % % % %     
% % % % % % % % % % %     for iz=1:nz
% % % % % % % % % % %          z=par.modz(iz);
% % % % % % % % % % %          
% % % % % % % % % % %          if z>300
% % % % % % % % % % %          
% % % % % % % % % % %          x1=xw-350;
% % % % % % % % % % %          x2=x1+550+z;
% % % % % % % % % % %             
% % % % % % % % % % %          for ix=1:nx
% % % % % % % % % % %             if par.modx(ix)>x1 && par.modx(ix)<x2
% % % % % % % % % % %             x=par.modx(ix);
% % % % % % % % % % %             ya=-0.103*z-200+(xw-2*x)*-.015;
% % % % % % % % % % %             y2=ya+100; y1=y2-350;   
% % % % % % % % % % %             
% % % % % % % % % % %             for iy=1:ny
% % % % % % % % % % %                 if par.mody(iy)>y1 && par.mody(iy)<y2
% % % % % % % % % % %                     mf(iy,ix,iz)=(amp); 
% % % % % % % % % % %                     alerta=1;
% % % % % % % % % % %                 end    
% % % % % % % % % % %             end
% % % % % % % % % % % 
% % % % % % % % % % %             end
% % % % % % % % % % % 
% % % % % % % % % % %          end
% % % % % % % % % % % 
% % % % % % % % % % %          end %(end if z <250)
% % % % % % % % % % %     end

 %%%   F%  Southern kaapvaal anomaly%%%
     xw=-100; ynw=400; 
    
    for iz=1:nz
         z=par.modz(iz);
         
         if  z>300 && z<600
         
         x1=xw-350;
         x2=x1+350+z;
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=-0.13*z-00+(xw-2*x)*.15;
            y2=ya+100; y1=y2-550;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=-amp*(-1); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
    end



% % %      
% % %      hold on
% % %   %%%%%%%%G!  %%%% Namaqua Natal and CapeFold Belts
     xw=150; ynw=-400;    
     for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw-700; %x1 position
         x2=x1+550; %-(z*0.0001); %x2 position
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=.269*z-450+(xw-x)*(0.2145);  %positioning y
            %ya=-100+(xw-x)*-3;
            y2=ya-00; y1=y2-600;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=(amp*0); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
     end
% %     
% % %      
% % % 
% % % %    %%H!   simple anomaly to cover the complicated signal 
% % % % 
    xw=150; ynw=600;    
     for iz=1:nz
         z=par.modz(iz);
         
         if z>=300 && z<600
         
         x1=xw+50; %x1 position
         x2=x1+550;%-z; %x2 position
            
         for ix=1:nx
            if par.modx(ix)>x1 && par.modx(ix)<x2
            x=par.modx(ix);
            ya=.0003*z+650+(xw-x)*(.000145);  %positioning y
            %and 0.0003 is the movement rate in y direction 
            y2=ya+200; y1=y2-500;   
            
            for iy=1:ny
                if par.mody(iy)>y1 && par.mody(iy)<y2
                    mf(iy,ix,iz)=(amp*0); 
                    alerta=1;
                end    
            end

            end

         end

         end %(end if z <250)
     end
     

    
    df = G*mf(:);
    %noise = randn(size(df))*0.2;% Gaussian noise, standard deviation 
    %df = df + noise; % Add the noise
    DR = speye(nmodel);
    dR = [df; zeros(size(LR,1),1)];
    GR = [G; DR*damp];
    dR = [dR; zeros(nmodel,1)];
    GR = [GR; LR*smooth];
    [m,flag,relres,iter,resvec] = lsqr( GR, dR, 1e-6, 100 );
    % ---- plot synthetic inversion ------------
    make_DS_map_synth_crt(mf,m,par,stax,stay); 
    upsample_model_syn_mf_crt(mf,par,1);
    toc
    upsample_model_syn_m_crt(m,par,1);
    toc
    
    
    
% % % %     %%%very nice anomaly in the middle of bushveld
% % % %     
% % % %     xw=-300; ynw=400; 
% % % %     
% % % %     for iz=1:nz
% % % %          z=par.modz(iz);
% % % %          
% % % %          if z<360
% % % %          
% % % %          x1=xw+350 ;
% % % %          x2=x1+350-z;
% % % %             
% % % %          for ix=1:nx
% % % %             if par.modx(ix)>x1 && par.modx(ix)<x2
% % % %             x=par.modx(ix);
% % % %             ya=-0.3*z-150+(xw-x)*(-.00015);
% % % %              y2=ya+650; y1=y2-300;   
% % % %             
% % % %             for iy=1:ny
% % % %                 if par.mody(iy)>y1 && par.mody(iy)<y2
% % % %                     mf(iy,ix,iz)=-amp; 
% % % %                     alerta=1;
% % % %                 end    
% % % %             end
% % % % 
% % % %             end
% % % % 
% % % %          end
% % % % 
% % % %          end %(end if z <250)
% % % %     end
% % % %     