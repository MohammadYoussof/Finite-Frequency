function dc = calc_dc(s_static,e_static,df,ray,stn)
disp('Mohammad, you know there is NO wonder that df = G*vm(:), is not?! ')
stn_nums = unique(ray.sta_num);

for ii = 1:length(stn_nums)
    indx = find(ray.sta_num==stn_nums(ii));
    df(indx) = df(indx)-s_static(ii);
end

orids = unique(ray.orid);
for ii = 1:length(unique(ray.orid))
    indx     = find(ray.orid==orids(ii));
    df(indx) = df(indx)-e_static(ii);
end

dc = df;