%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-2-3(yyyy-mm-dd)
% 贝塞尔多项式推导
%--------------------------------------------------------------------------
function [P, wr, absND] = fun_gaussian_polynomial(n, ep, varargin)
if nargin>2
    disp_en = varargin{1};
else
    disp_en = 0;
end
absND   = vpa(zeros(1, 2*n+1));
f = vpa(1);
for ii=1:n+1
    if ii==1
        f = 1;
    else
        f = vpa(f*(ii-1));
    end
    absND(2*ii-1) = vpa(1/f);
end
absNDr    = eval(absND);
absNDr(1) = absNDr(1)-(1+ep.^2);
rND       = roots(absNDr);
wr        = rND(real(rND)>0 & abs(imag(rND))<1e-5);
absND2    = absND;
absND2(3:4:end) = -absND2(3:4:end);% s=jw求解根
absP      = roots(fliplr(absND2));
P         = absP(real(absP)<0);
P         = P';
end