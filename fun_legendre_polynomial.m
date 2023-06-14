%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-2-3(yyyy-mm-dd)
% 勒让德多项式推导(Legendre Polynomial)
%--------------------------------------------------------------------------
function [P, wr, absND] = fun_legendre_polynomial(n, ep, varargin)
if nargin>2
    disp_en = varargin{1};
else
    disp_en = 0;
end
[Ln] = funGenLegendreCharacteristicPoly(n);
absND = Ln;
absND(1) = 1;
absNDr    = absND;
absNDr(1) = absNDr(1)-(1+ep.^2);
rND       = roots(absNDr);
wr        = rND(real(rND)>0 & abs(imag(rND))<1e-5);
absND2    = absND;
absND2(3:4:end) = -absND2(3:4:end);% s=jw求解根
absP      = roots(fliplr(absND2));
P         = absP(real(absP)<0);
P         = P';
end