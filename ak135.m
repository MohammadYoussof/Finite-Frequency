function varargout = ak135( varargin )

persistent discons;
persistent Re;
persistent VpCoeffs;
persistent VpdCoeffs;
persistent Vpd2Coeffs;
persistent VsCoeffs;
persistent VsdCoeffs;
persistent Vsd2Coeffs;

if ( isempty(VpCoeffs) )
  ak135_coeff;
end

errStr = nargchk( 1, 2, nargin );
if (~ isempty(errStr) )
  error( errStr );
end

errStr = nargoutchk(0,6,nargout);
if (~ isempty(errStr) )
  error( errStr );
end

if ( nargin == 1 && strcmp('discons',lower(varargin{1})) )
  varargout{1} = discons(1:end-1);
  return;
else
  r = varargin{1};
end

if ( any(r < 0) )
  error( 'ak135(): Radius must be positive.' );
end

if ( nargin == 2 )
  layer = varargin{2};
  if ( any(r > discons(layer) | r < discons(layer+1)) ) 
    error( 'ak135():  Out of layer.' );
  end
end

%if ( any(ismember(discons(2:end-1),r)) )
%  if ( nargin == 2 )
%    layer = varargin{2};
%  else
%    error( 'ak135(): Specify layer for discontinuities.' );
%  end
%end

nout = max([1 nargout]);
varargout = cell( 1, nout );
[varargout{:}] = deal( zeros(size(r))*NaN );


for d = 1:length(discons)-1
  
  indx = find( r <= discons(d) & r >= discons(d+1) );
  
  if ( isempty(indx) )
    continue;
  end
  
  if ( nargin == 2 && d ~= layer )
    continue;
  end
  
  x = r(indx);
  
  varargout{1}(indx) = polyval( VpCoeffs(d,:), x );  
  
  if ( nargout >= 2 )
    varargout{2}(indx) = polyval( VsCoeffs(d,:), x );
  end
  
  if ( nargout >= 3 )
    varargout{3}(indx) = polyval( VpdCoeffs(d,:), x );
  end
  
  if ( nargout >= 4 )
    varargout{4}(indx) = polyval( VsdCoeffs(d,:), x );
  end
  
  if ( nargout >= 5 )
    varargout{5}(indx) = polyval( Vpd2Coeffs(d,:), x );
  end
  
  if ( nargout == 6 )
    varargout{6}(indx) = polyval( Vsd2Coeffs(d,:), x );
  end
  
end

if ( any(isnan([varargout{:}])) )
  warning( 'ak135():  Only partial answer (filled with NaNs).' );
end


function out = akPoly( coeffs, r )
  
if ( size(r,2) == 1 )
  r = r';
end

exponents = ( size(coeffs,2)-1:-1:0 )';
rMatrix = repmat( r, length(coeffs), 1 );
exponentsMatrix = repmat( exponents, 1, length(r) );
rMatrix = rMatrix .^ exponentsMatrix;

out = coeffs * rMatrix;

