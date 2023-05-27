%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-11(yyyy-mm-dd)
% 多项式展开
% Fs = [a0, a1, a2, a3, ...]
%       a0 + a1*w + a2*w^2 + a3*w^3 + ...
%--------------------------------------------------------------------------
function Fs = funRecursionPoly(n, Zv)
    Fs(1) = vpa(1);
    nn = length(Zv);
    for m=1:nn
        zm = Zv(m);
        Fs(m+1) = Fs(m);
        for k = m :-1: 2
            Fs(k) = Fs(k-1) - Fs(k)*vpa(zm);
        end
        Fs(1) = -Fs(1)*vpa(zm);
    end
    if n>nn
        for ii=nn+2:n+1
            Fs(ii) = 0;
        end
    end