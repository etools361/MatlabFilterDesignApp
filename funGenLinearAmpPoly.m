%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-30(yyyy-mm-dd)
% 线性幅度滤波器系数设计
% 注意返回值低次在前，高次在后
%--------------------------------------------------------------------------
function [Ln] = funGenLinearAmpPoly(FilterOrder, IL)
if mod(FilterOrder, 2)
    Ln = [];
    warning('FilerOrder must be even');
    return;
end
% 使用chebfun函数计算m和b，使得在w=1处有单位增益和最高点
[m, b] = funLinearAmpGet_mb_chebfun(IL, FilterOrder);
if FilterOrder > 20
    dig = 100;
else
    dig = 16;
end
syms x
% 定义了直线幅度滤波器传递函数
f = 1/((1-IL)*sqrt(x)+b)^2;
a = sym(m);
coeffsx = sym(0);
for k = 0:FilterOrder
    kk = sym(k);
    % 使用在x=m处的麦克劳林展开得到系数
    coeffsx = coeffsx + vpa(subs(diff(f, x, kk) / factorial(k), x, a)*vpa((x-a)^kk, dig), dig);
end
polyx = coeffs(coeffsx);
Ln = vpa(zeros(1, length(polyx)*2-1), dig);
Ln(1:2:end) = polyx;
end



