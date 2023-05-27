%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-2-12(yyyy-mm-dd)
% 求H(jw)^2
%--------------------------------------------------------------------------
function [absH] = funCalcuHjw2(H)
n    = length(H);
Xpjw    = zeros(1, n);
Xnjw    = zeros(1, n);
Xpjw(1) = 1;
Xnjw(1) = 1;
for ii=2:n
    Xpjw(ii) = 1i^(ii-1);
    Xnjw(ii) = (-1i)^(ii-1);
end
HX1 = H.*Xnjw;
sHX1 = zeros(1, 2*n);
sHX1(n+1:end) = HX1;
HX2 = fliplr(H.*Xpjw);
for ii=1:2*n-1
    sHX1 = [sHX1(2:end),0];
    absH(ii) = sum(sHX1(1:n).*HX2);
end
% absH = conv((H.*Xpjw), (H.*Xnjw));
end