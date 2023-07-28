%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-6-30(yyyy-mm-dd)
% 线性幅度多项式推导
%--------------------------------------------------------------------------
function [P, wr, absND] = fun_linear_amp_polynomial_b(n, b, varargin)
if nargin>2
    disp_en = varargin{1};
else
    disp_en = 0;
end
[Ln] = funGenLinearAmpPoly_b(n, b);
if isempty(Ln)
    P = [];
    wr = 1;
    absND = [];
    return;
end
absND     = Ln;
wr        = 1;
absND2    = Ln;
absND2(3:4:end) = -absND2(3:4:end);% s=jw求解根
absP      = roots(fliplr(absND2));
P         = absP(real(absP)<0);
P         = P';
end