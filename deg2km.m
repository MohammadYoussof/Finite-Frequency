function km = deg2km( deg, r )

if ( nargin == 1 )
    r = 6371;
end

km = rad2km( deg2rad(deg), r );

