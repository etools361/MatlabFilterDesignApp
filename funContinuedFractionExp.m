%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-11(yyyy-mm-dd)
% 连分式展开
%--------------------------------------------------------------------------
function km = funContinuedFractionExp(n, Z, P, varargin)
% 辗转相除算法
km = [];
for ii=1:n
    ii_u = n+2-ii;
    if abs(Z(ii_u))<1e-12
        Z1    = [0,Z(1:ii_u-1),zeros(1,ii-1)];
    else
        Z1    = Z;
    end
    km(ii) = P(ii_u)/Z1(ii_u);
    P = P-Z1.*P(ii_u)/Z1(ii_u);
    T = P;
    P = Z;
    Z = T;
end
km = abs(km);
