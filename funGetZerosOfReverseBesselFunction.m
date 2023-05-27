%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% reverse bessel zeros
%--------------------------------------------------------------------------

x0 = 1;
a  = 2;
n  = 10;
epslion = 1e-10;
Delta = 1 + epslion;
z  = x0.*sqrt((n+a/2)*(n+a/2-1));

while Delta>epslion
    y = z;
    Omega = -1+(2-a)./z-(n+a/2)*(n+a/2-1)./z.^2;
    Tz = z-1./sqrt(Omega).*atan(sqrt(Omega).*theta_z_a./d_theta_z_a);
    z = Tz;
    Delta = abs(z-y)./abs(y);
end

