%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-02(yyyy-mm-dd)
% Chebyshev 滤波器综合，实现了低通原型参数的计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km] = funSynthesisInverseChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As)
    if isempty(Ap) || Ap<0
        Ap = 3;
        fprintf('Ap=%f dB\n', Ap);
    end
    if isempty(As) || As<0
        As = 70;
        fprintf('As=%f dB\n', As);
    end
    % Inverse Chebyshev参数计算
    if isempty(n) || n < 2
        % 由 As计算n
        n_min = acosh(sqrt((10^(-0.1*As)-1)/(10^(0.1*Ap)-1)))/acosh(fs/fp);
        n = ceil(n_min);
        fprintf('Order=%d\n', n);
    end
    [cellValueNetlist, km] = funEvenOrderParameter(n, Rs, Rl, Ap, As, fp);

function [cellValueNetlist, km] = funEvenOrderParameter(n, Rs, Rl, Ap, As, fp)
%     if Rs == Rl
%         Rl = Rl*(1+1e-6);
%     end
    if Rs>Rl
        t = sqrt(Rs/Rl);
    else
        t = sqrt(Rl/Rs);
    end
    % calcu Fs
    epsilon   = 1/sqrt(10^(0.1*As)-1);% 阻带衰减量
    epsilon2  = 1/sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
    phi2      = 1/n*asinh(1/epsilon);
    K         = cosh(1/n*acosh(1/(epsilon/epsilon2)));
    v2        = (n-1)*pi/(2*n);
    if ~mod(n, 2)
        K = sqrt((K).^2 + cos(v2).^2)./sin(v2);
    end
    if Rs == 0 || Rs == inf || Rl == 0 || Rl == inf
        B = 0;
    else
        B         = 1-((t^2-1)/(t^2+1))^2;
    end
    m         = (1-B)*(epsilon/epsilon2)^2;
    kk        = (1+2/m+2*sqrt(m+1)/m)^(1/2/n);
    h         = (sqrt(1+m)+sqrt(m))^(1/n);
%     kk = h;
    Zv        = zeros(1, n);
    Rv        = zeros(1, n);
    Tn        = zeros(1, n);
%     K       = (sqrt((K).^2 + cos(v2).^2)./sin(v2));
    for ii=1:n
        k  = ii;
        v  = (2*k-1)*pi/(2*n);
        Zv(ii)  = (-sinh(phi2).*sin(v) + 1i.*cosh(phi2).*cos(v)); % 8.52
        Rv(ii)  = -1i.*cos(v);%0.5*(h*exp(1i*(pi/2+v))+1/h*exp(1i*(pi/2-v)));%
        Tn(ii)  = 0.5*(kk*exp(1i*(pi/2+v))+1/kk*exp(1i*(pi/2-v)));%-1i.*cos(v);
    end
    if mod(n, 2)
        Zv1 = K./Zv;
        Rv1 = K./Rv;
        Tn1 = K./Tn;
    else
        Zv1 = K./(sqrt((Zv).^2 + cos(v2).^2)./sin(v2));
        Zv1 = -abs(real(Zv1))+1i.*imag(Zv1); % 可能存在问题
        Rv1 = K./(sqrt((Rv).^2 + cos(v2).^2)./sin(v2));
%         Rv1 = K./Rv;
        Rv1 = -abs(real(Rv1))+1i.*imag(Rv1); % 可能存在问题
%         Tn1 = K./Tn;
        Tn1 = K./(sqrt((Tn).^2 + cos(v2).^2)./sin(v2));
        Tn1 = -abs(real(Tn1))+1i.*imag(Tn1); % 可能存在问题
    end
    for jj=1:(mod(n+1,2)+1)
        iRv1 = abs(Rv1); % 找出最大极点并丢掉
        imax = max(iRv1)==iRv1;
        Rv1(imax) = [];
        if sum(imax) == 2
            break;
        end
    end
%     Rv1(abs((Rv1))>100*fp*2*pi) = [];% 排除掉无用点
%     Rv1     = imag(Rv1).*1i;
    Es      = funRecursionPoly(n, Zv1);
    Es      = abs((Es)); % 系数为正实数
    Ps      = funRecursionPoly(n, Rv1);
    Ps      = abs((Ps)); % 系数为正实数
    Fs      = funRecursionPoly(n, Tn1);
    Fs      = abs((Fs)); % 系数为正实数
%     Fs      = zeros(1, n+1);
%     Fs(n+1) = 1;
    if Rs == 0 || Rs == inf || Rl == 0 || Rl == inf
        % 一端接载
        Z  = Es;
        P  = Es;
        if mod(n, 2)
            Z(2:2:end) = 0;
            P(1:2:end) = 0;
        else
            Z(1:2:end) = 0;
            P(2:2:end) = 0;
        end
    else
        % 两端接载
        Z  = Es-Fs;
        P  = Es+Fs;
        % 系数归一化
    end
    [P, Z] = funPolyHighestOrderNorm(P, Z);
    % 零点移位法(适用于带零点的滤波器综合)
    % 求零极点
    if abs(real(Rv1))>1e-3
        for ii=1:length(Rv1)
            fprintf('Zeros%d=%0.3f+%0.3fi\n', ii, real(Rv1(ii)), imag(Rv1(ii)));
        end
    end
    nz   = ceil(length(Rv1)/2);
    ZAll = Rv1(1:nz);%1i.*imag(Rv1(imag(Rv1)<0));
    ZRm  = ZAll;
    nRmZ = length(ZRm);
    Z1P  = P;
    Z1Z  = Z;
    cn   = 1;cn2 = 1;
    % 对于偶数阶chebyshev II，需要先综合一个电感，后续和奇数阶相同
    cellValueNetlist = [];
    if ~mod(n, 2)
        Z1P  = Z;
        Z1Z  = P;
        Z1Ps   = funPolyMut_s(Z1P);
        km(cn) = Z1Z(end)./Z1Ps(end);
        cellValueNetlist{cn2} = {'L', 'S', real(km(cn))}; cn2 = cn2 + 1;
        ZTemp  = Z1Z-Z1Ps.*km(cn);
        Z1Z    = ZTemp;
        cn     = cn + 1;
        [Z1P, Z1Z] = funPolyHighestOrderNorm(Z1P, Z1Z);
    end
    % PI型
    for ii = 1:nz
        Z1Zs = funPolyMut_s(Z1Z);
        % 找电容最小的零点
        C1Temp = [];
        for jj=1:nRmZ
            C1Temp(jj)   = funGetPolyValue(Z1P, ZRm(jj))/funGetPolyValue(Z1Zs, ZRm(jj));
        end
        [iC1] = find(C1Temp>0);
        if ~isempty(iC1)
            [a, b] = min(C1Temp(iC1));
            iMin   = iC1(b);
            C1     = a;
        else
            iMin   = 1;
            C1     = C1Temp(iMin);
        end
        ZRm0   = ZRm(iMin); % 最终选择的所要去除的0点
        ZRm(iMin) = [];
        nRmZ      = nRmZ - 1;
        if nRmZ < 0
            nRmZ = 0;
        end
        %---
        km(cn) = C1;cn = cn + 1;
        cellValueNetlist{cn2} = {'C', 'P', real(C1)};cn2 = cn2 + 1;
        Y3P    = Z1Z;
        Y3Z    = Z1P-C1.*Z1Zs;
        [Y3P, Y3Z] = funPolyHighestOrderNorm(Y3P, Y3Z);
        % 求wi=ZAll(nz+1-ii)点的留数
        Po = ZRm0;
        Po = conj(Po)*Po;
        Y3ZPr = funPolyDivPole(Y3Z, Po);
        Ki = funGetPolyValue(Y3P, ZRm0)/funGetPolyValue(funPolyMut_s(Y3ZPr), ZRm0);
        % 求器件值
        C2 = 1/Ki;
        L2 = 1/(C2*Po);
        if ~mod(n, 2)
            km(cn) = C2;cn = cn + 1;
            km(cn) = L2;cn = cn + 1;
        else
            km(cn) = L2;cn = cn + 1;
            km(cn) = C2;cn = cn + 1;
        end
        cellValueNetlist{cn2} = {'LCP', 'S', real(L2), real(C2)};cn2 = cn2 + 1;
        % 求去除掉L2和C2后的阻抗
        Z4P = Y3ZPr;
        Z4Z = funPolyDivPole(Y3P-Ki*funPolyMut_s(Y3ZPr), Po);
        [Z1P, Z1Z] = funPolyHighestOrderNorm(Z4P, Z4Z);
    end
    if nz == 0
        Z4Z = Z1Z;
    end
    C3  = 1/Z4Z(1);
    km(cn) = C3;
    cellValueNetlist{cn2} = {'C', 'P', real(C3)};cn2 = cn2 + 1;
    if ~mod(n, 2) && Rl==Rs
        km  = fliplr(km);
    end
    cellValueNetlist = fliplr(cellValueNetlist);
    km  = real(km);
    
    % 辗转相除算法(仅仅适用于全极点滤波器综合)
    % km = funContinuedFractionExp(n, Z, P);


