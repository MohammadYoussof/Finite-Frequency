function yesno = isEven( x )

% $Id: isEven.m,v 1.1 2006/01/06 00:27:32 jason Exp jason $

    yesno = mod( x, 2 ) == 0;
    