%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-13(yyyy-mm-dd)
% Legendre 滤波器综合，实现了低通原型参数计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km, R1] = funSynthesisLinearAmpFilter(n, Rs, Rl, fp, fs, Ap, As)
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
    [cellValueNetlist, km, R1] = funEvenOrderParameter(n, Rs, Rl, Ap);

function [cellValueNetlist, km, R1] = funEvenOrderParameter(n, Rs, Rl, Ap)
    if n>10
        digits(40)
    else
        digits(32)
    end
    if Rs == Rl
        Rl = Rl*(1+1e-12);
    end
    if Rs>Rl
        t = sqrt(Rs/Rl);
    else
        t = sqrt(Rl/Rs);
    end
%     epsilon   = sqrt(10^(0.1*Ap)-1);
    % calcu Fs
    IL   = 10^(-Ap/20);
    [P, w1, absND] = fun_linear_amp_polynomial(n, IL);
    epsilon_x = 1/sqrt(absND(1));
    RL = epsilon_x/(1-epsilon_x)*Rs;
    ND     = (funRecursionPoly(n, P));
    NDflip = ND;
    Fs = abs(real(NDflip));
    Fe  = Fs;
    Fe(1:2:end)  = 0;
    Fo  = Fs;
    Fo(2:2:end)  = 0;
    RX = RL/(Rs+RL);
    Z = Fe;
    P = Fo.*RX;
    % 辗转相除算法
    km = funContinuedFractionExp(n, Z, P);
    km = fliplr(km);
    R1 = RL;
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