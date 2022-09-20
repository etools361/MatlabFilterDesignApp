%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-11(yyyy-mm-dd)
% 求多项式(P(nP)*s^(nP-1)+P(nP-1)*s^(nP-2)+...+P(2)*s^1+P(1)*s^0)和一个极点(s^2-Po)之比为:
% Pr(nP-2)*s^(nP-3)+Pr(nP-3)*s^(nP-4)+...+Pr(2)*s^1+Pr(1)*s^0
% 也即[P(nP), P(nP-1), ..., P(2), P(1)] = [Pr(nP-2),Pr(nP-3),...,Pr(2),Pr(1),0,0]  + [0,0,Po*Pr(nP-2),Po*Pr(nP-3),...,Po*Pr(2),Po*Pr(1)]
% 有：
% P(nP)   = Pr(nP-2);                -->  Pr(nP-2) = P(nP)
% P(nP-1) = Pr(nP-3);                -->  Pr(nP-3) = P(nP-1)
% P(nP-2) = Pr(nP-4) + Po*Pr(nP-2);  -->  Pr(nP-4) = P(nP-2) - Po*Pr(nP-2)
% P(nP-3) = Pr(nP-5) + Po*Pr(nP-3);  -->  Pr(nP-5) = P(nP-3) - Po*Pr(nP-3)
% ...
% P(4)    = Pr(2)    + Po*Pr(4);     -->  Pr(2)    = P(4)    - Po*Pr(4)
% P(3)    = Pr(1)    + Po*Pr(3);     -->  Pr(1)    = P(3)    - Po*Pr(3)
% P(2)    =            Po*Pr(2);     -->  用于验证
% P(1)    =            Po*Pr(1);     -->  用于验证
%--------------------------------------------------------------------------
function [Pr] = funPolyDivPole(P, Po)
nP = length(P);
Pr = zeros(1, nP);
for ii=1:nP-2
    if ii == 1 || ii == 2
        Pr(nP-2+1-ii) = P(nP+1-ii);
    else
        Pr(nP-2+1-ii) = P(nP+1-ii)-Po*Pr(nP+1-ii);
    end
end
% 验证
MinErr = 1e-4;
if abs(P(2)-Po*Pr(2))<MinErr && abs(P(1)-Po*Pr(1))<MinErr
%     fprintf('OK\n');
else
    fprintf('Error! Can not converge(%e, %e)\n', abs(P(2)-Po*Pr(2)), abs(P(1)-Po*Pr(1)));
%     Pr = [];
end
end