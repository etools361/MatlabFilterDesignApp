%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-02-18(yyyy-mm-dd)
% Gaussian 滤波器综合，实现了低通原型参数计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km] = funSynthesisGaussianFilter2(n, Rs, Rl, fp, fs, Ap, As)
    if isempty(Ap) || Ap<0
        Ap = 3;
        fprintf('Ap=%f dB\n', Ap);
    end
    % Chebyshev参数计算
    if isempty(n) || n < 2
        % 由 As计算n
        n_min = 1/acosh(fs/fp)*acosh(sqrt((10^(-0.1*As)-1)/(10^(0.1*Ap)-1)));
        n = ceil(n_min);
        fprintf('Order=%d\n', n);
    end
    [cellValueNetlist, km] = funEvenOrderParameter(n, Rs, Rl, Ap);

function [cellValueNetlist, km] = funEvenOrderParameter(n, Rs, Rl, Ap)
    if Rs == Rl
        Rl = Rl*(1+1e-6);
    end
    if Rs>Rl
        t = sqrt(Rs/Rl);
    else
        t = sqrt(Rl/Rs);
    end
    eta = ((t+1/t)/2)^2;
    epsilon   = sqrt(10^(0.1*Ap)-1);
    % calcu Fs
    [P, w1, absND] = fun_gaussian_polynomial2(n, epsilon);
    [absH] = absND;
    ND     = fliplr(poly(P));
    if isinf(eta)
        K2    = absH;
        K2(1) = 0;
    else
        K2    = eta.*absH./absH(1);
        K2(1) = K2(1)-1;
    end
    K2(abs(K2)<1e-10)=0;
    m = length(K2);
    KK = zeros(1, m);
    for ii=1:m
        KK(ii) = 1i.^(ii-1);
    end
    K2 = K2.*KK;
    rK2 = roots(fliplr(K2));
    absrK2 = abs(rK2);
    [a, b]=sort(absrK2);
    nZeros = 0;
    for ii=1:m
        if a(ii)>1e-5
            break;
        end
        nZeros = nZeros + 1;
    end
    rootsSel = rK2(b(ii:end));
    rootsSel = rootsSel((abs(real(rootsSel))>1e-10&real(rootsSel)<0) | (abs(real(rootsSel))<1e-10));
    nZeroAdd = floor(nZeros/2);
    for ii=1:nZeroAdd
        rootsSel(end+1) = 0;
        if length(rootsSel)>=n
            break;
        end
    end
    Fsflip = poly(rootsSel.*w1);
    rND    = roots(fliplr(ND));
    NDflip = poly(rND.*w1);
    Fs = fliplr(NDflip);
    Es = fliplr(Fsflip);
    if Rs == 0 || Rs == inf || Rl == 0 || Rl == inf
        % 一端接载
        Z  = Fs;
        P  = Fs;
        if mod(n, 2)
            Z(2:2:end) = 0;
            P(1:2:end) = 0;
        else
            Z(1:2:end) = 0;
            P(2:2:end) = 0;
        end
    else
        Ee  = Es;
        Ee(1:2:end)  = 0;
        Eo  = Es;
        Eo(2:2:end)  = 0;
        Fe  = Fs;
        Fe(1:2:end)  = 0;
        Fo  = Fs;
        Fo(2:2:end)  = 0;
        if mod(n, 2)
            Z   = Eo-Fo;
            P   = Ee+Fe;
        else
            P   = Eo+Fo;
            Z   = Ee-Fe;
        end
        % 两端接载
%         Z  = Fs-Es;
%         P  = Fs+Es;
    end
    % 辗转相除算法
    km = funContinuedFractionExp(n, Z, P);
    km = fliplr(km);
    cellValueNetlist = [];
    for ii=1:n
        if mod(ii, 2)
            Type = 'C';
            SP   = 'P';
        else
            Type = 'L';
            SP   = 'S';
        end
        cellValueNetlist{ii} = {Type, SP, km(ii)};
    end