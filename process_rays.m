%function [data, raw_nrays, reduced_nrays] = process_delays(data)
clear all
close all

%=================== input rays ========================================
fid = fopen('n3.dat', 'r');
A = [textscan(fid, '%f %f %f %s %d %f %f %f %f %d'),1]; fclose(fid);
ray.pd    = A{1};
ray.p     = ray.pd .*(360/(2*pi*6371));
ray.baz   = A{2};
d         = A{3};
ray.d     = d;
ray.sta   = A{4};
ray.orid  = A{5};
ray.cf    = A{6};
ray.olat  = A{7};
ray.olon  = A{8};
ray.odep  = A{9};
ray.chan  = A{10};
sta_ray   = ray.sta;
ray.nrays = length(ray.p);

%=================== input stations ==================================== 
sta_names = unique( ray.sta );
fid = fopen('new_crust_model.txt','r');
B = [textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f'),1]; fclose(fid); 
stns.sta = B{1};  stns.lat = B{2};  stns.lon = B{3};  stns.elv = B{4};
stns.z   = B{5};  stns.n   = B{6};  stns.ne  = B{7};  stns.e   = B{8};
stns.se  = B{9};  stns.s   = B{10}; stns.sw  = B{11}; stns.w   = B{12}; 
stns.nw  = B{13}; 
stns.sta = unique(stns.sta);

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
% ------- associate data ------------------------------------------------
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

%---------- site corrections ----------------------------------------------
data = crust_corr(data);



%=========================================================================

% set bin increments
baz_dk = 0:8:360;  ss_dk = [4 5.5 6.4 7.5 10];

ray = data.ray; stn = data.stn; 

ray.ss = ray.pd;
ray.m1_diff = zeros(length(ray.ss),1); ray.flag = ones(length(ray.ss),1);
ray.m2_diff = zeros(length(ray.ss),1); ray.f_band = zeros(length(ray.ss),1);
ray.d2 = ray.d - ray.corr1;

band1 = find(ray.cf==1);
ray.f_band(band1) = 1;
band2 = find(ray.cf==0.5);
ray.f_band(band2) = 2;
band3 = find(ray.cf==0.3);
ray.f_band(band3) = 3;



% Processing rays by baz, ss
h = waitbar( 0, 'Processing Rays for reducing the major trends...' );
for cnt=1:length(stn.num)
    waitbar( (cnt-1)/length(stn.num), h );
    for i_fband = 1:max(ray.f_band)
    ind1 = find(ray.sta_num==stn.num(cnt));
        for baz_ind = 1:(length(baz_dk)-1)
            for ss_ind = 1:(length(ss_dk)-1)
                ind2a = find(ray.baz >= baz_dk(baz_ind));
                ind2b = intersect(ind1,ind2a);
                ind3  = find(ray.baz < baz_dk(baz_ind+1));
                ind23 = intersect(ind2b,ind3);
                ind4a = find(ray.ss >= ss_dk(ss_ind));
                ind4b = intersect(ind1,ind4a);
                ind5  = find(ray.ss < ss_dk(ss_ind+1));
                ind45 = intersect(ind4b,ind5);
                bin_indx = intersect(ind23,ind45);
                bin_n = length(bin_indx);
            
                if  ( bin_n>=5 && bin_n<7 )
                    bin_m1 = mean(ray.d2(bin_indx));
                    ray.m1_diff(bin_indx) = abs(ray.d2(bin_indx)-bin_m1);
                    limit = max(ray.m1_diff(bin_indx));
                    for jj = 1:length(bin_indx)
                       if ray.m1_diff(bin_indx(jj))>=limit
                           ray.flag(bin_indx(jj))=0;
                       else
                           ray.flag(bin_indx(jj))=1;
                       end
                    end
                elseif bin_n>6 
                    bin_m1 = mean(ray.d2(bin_indx));
                    ray.m1_diff(bin_indx) = abs(ray.d2(bin_indx)-bin_m1);
                    limit = max(ray.m1_diff(bin_indx));
                    bin_total = 0;
                    for jj = 1:length(bin_indx)
                       if ray.m1_diff(bin_indx(jj))>=limit
                           ray.flag(bin_indx(jj))=0;
                       else
                           ray.flag(bin_indx(jj))=1;
                           bin_total = bin_total + ray.d2(bin_indx(jj));
                       end
                    end
                    bin_m2 = bin_total/sum(ray.flag(bin_indx));
                    ray.m2_diff(bin_indx) = abs(ray.d2(bin_indx)-bin_m2);
                    m2_diff = sort(ray.m2_diff(bin_indx));
                    limit = m2_diff(5);
                    for kk = 1:length(bin_indx)
                       if ray.m1_diff(bin_indx(kk))>limit
                           ray.flag(bin_indx(kk))=0;
                       else
                           ray.flag(bin_indx(kk))=1;
                       end
                    end
                elseif bin_n < 3
                    ray.flag(bin_indx)=1;
                end 
            end
        end
    end
end      
close( h );

raw_nrays = length(ray.ss)
reduced_nrays = sum(ray.flag)

use_rays = find(ray.flag==1);
ray2.sta=ray.sta(use_rays);          ray2.d=ray.d(use_rays);              ray2.baz = ray.baz(use_rays); 
ray2.pd = ray.pd(use_rays);          ray2.cf=ray.cf(use_rays);    
ray2.orid=ray.orid(use_rays);       
ray2.olat=ray.olat(use_rays);        ray2.olon=ray.olon(use_rays);        ray2.odep=ray.odep(use_rays);
ray2.chan=ray.chan(use_rays);    

%chan =1;

fid = fopen('TDelay_Reduced_Evs.txt','w');
for qq=1:length(ray2.d)
       sta = char(ray2.sta(qq));
       fprintf(fid,' %f   %f   %f   %s   %d   %f   %f   %f   %f   %d \n'...
           ,ray2.pd(qq),ray2.baz(qq),ray2.d(qq),sta...
           ,ray2.orid(qq),ray2.cf(qq),ray2.olat(qq),ray2.olon(qq),ray2.odep(qq),ray2.chan(qq));
        
end
fclose(fid);  

%data.ray = ray;
%data.stn = stn;


