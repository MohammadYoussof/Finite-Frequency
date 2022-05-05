function K = Kval(par,a,rp,theta,kp,p,baz,cf,Vz,chan)

%disp('wavelength from x-corr center frequency ')
%disp('1st Fresnel, calibrated from BD kernels')

modz = par.modz; z_ind = kp;  zp = modz(kp);
if zp == modz(1)
    dz = modz(2)-modz(1);
else
    dz = 0.5*( (modz(z_ind+1)-modz(z_ind))+(modz(z_ind)-modz(z_ind-1)) );
end
modx = par.modx; mody = par.mody; x_ind = find(modx == a(1)); y_ind = find(mody == a(2));
if ( x_ind>3 && x_ind<(length(modx)-2) && y_ind>3 && y_ind<(length(mody)-2) )
    dh = 0.25*( (modx(x_ind+1)-modx(x_ind))+(modx(x_ind)-modx(x_ind-1)) ...
                +(mody(y_ind-1)-mody(y_ind))+(mody(y_ind)-mody(y_ind+1)) );
elseif ( x_ind<=3 || x_ind>=(length(modx)-2) || y_ind<=3 || y_ind>=(length(mody)-2) )
    dh = modx(2)-modx(1);
end

r_scal = (6371-zp)/6371;
v   = Vz(round(zp));
wl  = v/cf;                                            % ~wavelength from x-corr center frequency  
F1r = sqrt(0.235*wl*(zp/cos(asin(p*Vz(275)))));         % ~1st Fresnel, calibrated from BD kernels
% v
% disp('wavelength')
% wl
% disp('Frensel')
% F1r

K = zeros(size(rp));
rfe1  = (90-baz)*(pi/180);
xb = cos(rfe1);
yb = sin(rfe1);
rfe = atan2(yb,xb);
phi  = -theta + rfe;
inc = asin(p*v);
if inc>1.4
    inc=1.4;
end
angwt = zeros(size(rp));
qq = 1:length(rp);
angwt(qq) = 1 + 0.6*(cos(phi(qq)).^2)*( (1/sin( (pi/2)-inc ))-1) ;
R_n1 = F1r*angwt;
phi = phi';
    if chan==0
    R_n1   = (0.85*dh)*angwt*(1+1.5*(zp/1400));
    use = rp./R_n1;
    ind = find(use<1);
    wt = ones(size(use(ind))) - use(ind).^2;
    wt = wt./sum(wt);  
    wt = ((-dz/cos(inc))/v)*r_scal*wt; 
    K(ind(:)) = wt(:); 
else  
    use = rp./R_n1;
    ind = find(use<1);
    if length(ind)<2
        R_n1 = R_n1 + (0.85*dh-F1r)*ones(size(R_n1));
        use = rp./R_n1;
        ind = find(use<1);
    end
    wt1 = (1.82-1.155*(rp(ind)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)+0.25*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt2 = (1.82-1.155*((rp(ind)+0.075*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)+0.075*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt3 = (1.82-1.155*((rp(ind)-0.075*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)-0.075*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt4 = (1.82-1.155*((rp(ind)+0.15*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)+0.15*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt5 = (1.82-1.155*((rp(ind)-0.15*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)-0.15*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt6 = (1.82-1.155*((rp(ind)+0.25*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)+0.25*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt7 = (1.82-1.155*((rp(ind)-0.25*dh)./R_n1(ind))).*(1 + 0.15*(cos(0.5*phi(ind)).^2)).*...
            sin(pi*(((rp(ind)-0.25*dh+0.3*R_n1(ind))./(1.25*R_n1(ind))).^2));
    wt = (1/7)*( wt1+wt2+wt3+wt4+wt5+wt6+wt7 );
    wt = wt./sum(abs(wt));
    wt = ((-dz/cos(inc))/v)*r_scal*wt;
    %disp('equatoin_6=Kf1(R_N,D_R,ω,∆)');
    %m% save wt wt
    K(ind(:)) = wt(:);
    %m% save K K
% % % % %     disp('K')
% % % % % K

fid = fopen('Kernels_DATA','a');
fprintf(fid , ' %f  %f  %f  %f  %f  %f  %f  %s \n', v, wl, F1r, mean(R_n1), mean(rfe), mean(phi), mean(inc), '*' )
% fprintf(fid , ' %f  %f  %f  %f  %f  %f  %f  %f  %f  %f \n', v, wl, F1r, rfe1, rfe, phi, inc, R_n1, wt, K )
fclose(fid);
    end

end
