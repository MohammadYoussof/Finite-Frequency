Re = 6371;
discons = [6371 6351 6336 6251 6161 5961 5711 5611 3631 3479.5 1217.5 0]';

VpCoeffs = [...
  0.0       0.0         0.0       5.8;          %1
  0.0       0.0         0.0       6.5;          %2
  0.0       0.0         -.74953   8.78541;      %3
  0.0       0.0       -17.69722   25.41389;     %4
  0.0       0.0       -23.25415   30.78835;     %5
  0.0       0.0       -21.40656   29.38906;     %6
  0.0       0.0       -16.93412   25.97074;     %7
-26.6083   51.9932    -41.1538    25.1481;      %8
  0.0       0.0        -1.47089   14.4881;     %9
  0.0     -13.67046     3.75665   10.03024;     %10
  0.0      -4.09689     0.0       11.27214 ];   %11

VpCoeffs(:,1) = VpCoeffs(:,1)./(Re^3);
VpCoeffs(:,2) = VpCoeffs(:,2)./(Re^2);
VpCoeffs(:,3) = VpCoeffs(:,3)./Re;

VsCoeffs = [...
0            0       0           3.46;
0            0       0           3.85;
0            0      -2.248585    6.716231;
0            0       -1.27420    5.75020;
0            0      -11.08552   15.23852;
0            0      -13.50652   17.71792;
0            0      -16.53147   20.77960;
-14.1080    27.8988 -7.634931    12.9303;
0            0      -1.58205     8.15015;
0            0       0           0;
0           -3.45241 0           3.50440];

VsCoeffs(:,1) = VsCoeffs(:,1)./(Re^3);
VsCoeffs(:,2) = VsCoeffs(:,2)./(Re^2);
VsCoeffs(:,3) = VsCoeffs(:,3)./Re;

dCoeffs = repmat([3 2 1],size(VpCoeffs,1),1);
VpdCoeffs = VpCoeffs(:,1:3) .* dCoeffs;
VsdCoeffs = VsCoeffs(:,1:3) .* dCoeffs;

d2Coeffs = repmat([2 1],size(VpCoeffs,1),1);
Vpd2Coeffs = VpdCoeffs(:,1:2) .* d2Coeffs;
Vsd2Coeffs = VsdCoeffs(:,1:2) .* d2Coeffs;

clear dCoeffs d2Coeffs Re;
