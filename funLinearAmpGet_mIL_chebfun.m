%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-22(yyyy-mm-dd)
% 计算线性幅度滤波器归一化最佳值
% 由斜率和阶数计算出m和IL值
% n:滤波器阶数，必须是大于0的偶数
% 使用chebfun函数实现
%--------------------------------------------------------------------------
function [m, IL] = funLinearAmpGet_mIL_chebfun(b, n)
if mod(n, 2) || n < 2
    warning('n must be a even number');
end
mr = [0.3,0.95];
ILr = [0.1,0.95];
m  = chebfun2(@(m,IL) m, [mr, ILr]);
IL  = chebfun2(@(m,IL) IL, [mr, ILr]);
Lr = funLinearAmpEq(m, b, IL);
Lr0  = Lr;
X1   = Lr0-1;
X2   = 0;
for ii=1:n
    dfn = diff(Lr, ii, 2);
    X1 = X1 +    dfn./factorial(ii).*(1-m).^ii;
    X2 = X2 + ii.*dfn./factorial(ii).*(1-m).^(ii-1);
end
rx1 = roots(X1);
rx2 = roots(X2);
[~, ix1] = size(rx1);
[~, ix2] = size(rx2);
x00 = linspace(-1, 1, 1000);
y00 = linspace(-1, 1, 1000);
roots_ret = [];
iroot     = 0;
for cx1 = 1:ix1
    x1 = real(rx1(:, cx1));
    y1 = imag(rx1(:, cx1));
    for cx2 = 1:ix2
        x2 = real(rx2(:, cx2));
        y2 = imag(rx2(:, cx2));
        [x0,y0,iout,jout] = intersections(x1(x00),y1(y00),x2(x00),y2(y00));
        if ~isempty(x0)
            iroot = iroot + 1;
            roots_ret(iroot, :) = [x0(1), y0(1)];
        end
    end
end
if iroot==1
    m  = roots_ret(1);
    IL = roots_ret(2);
else
    warning('can not find the roots');
    m = 0.5;
    IL = b;
end
% ro = roots(X1,X2);
% m = ro(1);
% b = ro(2);
% fprintf('b:%f, m:%f\n', b, m);
end