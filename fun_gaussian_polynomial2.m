%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-2-3(yyyy-mm-dd)
% 多项式推导
%--------------------------------------------------------------------------
function [P, wr, absND] = fun_gaussian_polynomial2(n, ep, varargin)
if nargin>2
    disp_en = varargin{1};
else
    disp_en = 0;
end
absND   = zeros(1, 2*n+1);
% for ii=1:n+1
%     absND(2*ii-1) = 1/factorial(ii-1);
% end
ND(1) = 1;
for ii=1:n
    ND(ii+1) = 1/factorial(ii);
end
ND = ND./ND(end);
absND = funCalcuHjw2(ND);

absNDr    = absND./absND(1);
absNDr(1) = absNDr(1)-(1+ep.^2);
rND       = roots(absNDr);
wr        = rND(real(rND)>0 & abs(imag(rND))<1e-5);
absND2    = absND;
absND2(3:4:end) = -absND2(3:4:end);% s=jw求解根
absP      = roots(fliplr(absND2));
P         = absP(real(absP)<0);
P         = P';
end