function varargout = ray_trace( v1, v2, z1, z2, p, slat, slon, baz)


spherical_earth = true;

if ( nargout == 3 )
	calc_time = true;
else
	calc_time = false;
end

v1 = v1(:);
v2 = v2(:);
z1 = z1(:);
z2 = z2(:);

if ( spherical_earth )
	r1 = 6371-z1;
	r2 = 6371-z2;
	z1 = 6371*log(6371./r1);
	z2 = 6371*log(6371./r2);
	v1 = (6371./r1).*v1;
	v2 = (6371./r2).*v2;
end

u1 = 1./v1;
u2 = 1./v2;

dv = v2-v1;
dz = z2-z1;
b = dv ./ dz;
	
const_indx = b == 0;

% degenerate case
if ( p == 0 )

	X = zeros( length(v1)+1, 1 );
	if ( calc_time )
		T = zeros( 1, length(v1) );
		T(const_indx) = 1./v2(const_indx) .* dz;
		T(~const_indx) = (1./dv(~const_indx))./(2*dz(~const_indx)) + ...
			u1(~const_indx).*dz(~const_indx);
		T = [0; cumsum(T)];
	end
	
else

	X = zeros( length(v1), 1 );
	X(const_indx) = constv_dist( u1(const_indx), dz(const_indx), p );
	X(~const_indx) = gradv_dist( b(~const_indx), u1(~const_indx), ...
		u2(~const_indx), p );

	if ( calc_time )
		T = zeros( 1, length(u1) );
		T(const_indx) = constv_time( u1(const_indx), dz(const_indx), p );
		T(~const_indx) = gradv_time( b(~const_indx), u1(~const_indx), ...
			u2(~const_indx), p, X(~const_indx) );
		T = [0; cumsum(T)];
		indx = ~isreal(T);
		T(indx) = NaN;
	end

	X = [0; cumsum(X)];
	indx = find( imag(X) );
	X(indx) = NaN;
	X = km2deg(X);  %
end

if ( nargout == 1 )
	varargout{1} = X;
end

if ( nargout == 2 || nargout == 3 )
	if ( spherical_earth )
		[lat,lon] = reckon( slat, slon, X, baz );
		varargout{1} = lat;
		varargout{2} = lon;
	else
		lon = cos( az2dfe(deg2rad(baz)) ) .* X;
		lat = sin( az2dfe(deg2rad(baz)) ) .* X;
		varargout{1} = lat + slat;
		varargout{2} = lon + slon;
	end
end

if ( nargout == 3 )
	varargout{3} = T;
end


%
function t = constv_time( u, dz, p )

t = u.^2 .* dz ./ eta(u,p);

%
function x = constv_dist( u, dz, p )

x = (p * dz) ./ eta(u,p);
%i = asin(p/u);
%x = tan(i) * dz;

%
function t = gradv_time( b, u1, u2, p, x )
eta1 = eta( u1, p );
eta2 = eta( u2, p );
t1 = 1./b .* ( log((u1+eta1)/p) - eta1./u1 );
t2 = 1./b .* ( log((u2+eta2)/p) - eta2./u2 );
t = (t1 - t2) + p*x;

%
function x = gradv_dist( b, u1, u2, p )

x1 = eta(u1,p) ./ (b.*u1*p);
x2 = eta(u2,p) ./ (b.*u2*p);

x = x1 - x2;

%
function val = eta( u, p )

val = sqrt( u.^2 - p.^2 );

