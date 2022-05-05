function s = ray_length( rayxyz )
disp('Ray Length calculations')
s = [0; cumsum(sqrt( sum((diff(rayxyz)).^2,2) ))];
save s s