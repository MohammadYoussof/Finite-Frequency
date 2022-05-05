function notify( text, par, verbose_level )

if ( par.verbose >= verbose_level )
    fprintf( [text '\n'] );
end
