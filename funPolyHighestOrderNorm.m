%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-11(yyyy-mm-dd)
% 最高次项归一化
%--------------------------------------------------------------------------
function [P, Z] = funPolyHighestOrderNorm(P, Z)
% 最高次项系数归一化
nP = length(P);
BreakEn = 0;
for ii=1:nP
    if (Z(nP+1-ii)==1 && P(nP+1-ii)==0) || (Z(nP+1-ii)==0 && P(nP+1-ii)==1) || (Z(nP+1-ii)==1 && P(nP+1-ii)==1)
        break;
    end
    if Z(nP+1-ii)~=1 && Z(nP+1-ii)~=0
        KH = Z(nP+1-ii);
        BreakEn = 1;
    elseif P(nP+1-ii)~=1 && P(nP+1-ii)~=0
        KH = P(nP+1-ii);
        BreakEn = 1;
    else
        continue;
    end
    P = P./KH;
    Z = Z./KH;
    if BreakEn
        break;
    end
end
end