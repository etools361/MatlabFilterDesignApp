%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-11(yyyy-mm-dd)
% 最高次项归一化
%--------------------------------------------------------------------------
function [P, Z] = funPolyHighestOrderNorm(P, Z)
% 最高次项系数归一化
nP = length(P);
BreakEn = 0;
Err = 1e-9;
for ii=1:nP
    if (abs(Z(nP+1-ii)-1)<Err && abs(P(nP+1-ii))<Err) || (abs(Z(nP+1-ii))<Err && abs(P(nP+1-ii)-1)<Err) || (abs(Z(nP+1-ii)-1)<Err && abs(P(nP+1-ii)-1)<Err)
        break;
    end
    if abs(Z(nP+1-ii)-1)>Err && abs(Z(nP+1-ii))>Err
        KH = Z(nP+1-ii);
        BreakEn = 1;
    elseif abs(P(nP+1-ii)-1)>Err && abs(P(nP+1-ii))>Err
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