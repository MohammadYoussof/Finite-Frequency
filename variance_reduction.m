function vr = variance_reduction( d, df )

% vr = variance_reduction( d, df )
%
% d = measured residuals
% df = residuals from forward problem (df = G*m);

vr = 100*(1 - (var( d - df ) / var( d )));