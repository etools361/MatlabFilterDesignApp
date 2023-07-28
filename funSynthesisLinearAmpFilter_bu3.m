%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-13(yyyy-mm-dd)
% LinearAmplitude滤波器综合，实现了低通原型参数计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km, Rs] = funSynthesisLinearAmpFilter(n, Rs, Rl, fp, fs, Ap, As)
    if isempty(Ap) || Ap<0
        Ap = 3;
        fprintf('Ap=%f dB\n', Ap);
    end
    % 参数计算
    if isempty(n) || n < 2
        % 由 As计算n
        n_min = 1/acosh(fs/fp)*acosh(sqrt((10^(-0.1*As)-1)/(10^(0.1*Ap)-1)));
        n = ceil(n_min);
        fprintf('Order=%d\n', n);
    end
    [cellValueNetlist, km, Rs] = funEvenOrderParameter(n, Rs, Rl, Ap);

function [cellValueNetlist, km, Rs] = funEvenOrderParameter(n, Rs, Rl, Ap)
    if n>10
        digits(40)
    else
        digits(32)
    end
    eta = inf;
    if Rs == Rl
        Rl = Rl*(1+1e-12);
    end
%     epsilon   = sqrt(10^(0.1*Ap)-1);
    % calcu Fs
    IL   = 10^(-Ap/20);
    [P, w1, absND] = fun_linear_amp_polynomial(n, IL);
    f0 = 1/sqrt(absND(1));
%     V0 = 1;
    if Rl == inf % P
        RL0 = f0/(1-f0)*Rs; 
        SP0 = 'P';
        Pos = 'Head';
    elseif Rs == inf % P
        RL0 = 1/(1/f0-1/Rl);
        SP0 = 'P';
        Pos = 'Tail';
    elseif Rs == 0 % S
        RL0 = Rl*(1/f0-1);
        SP0 = 'S';
        Pos = 'Tail';
    elseif Rl == 0 % S
        RL0 = 1/f0-Rs;
        SP0 = 'S';
        Pos = 'Head';
    else
        R1 = (sqrt(5)-1)/2*Rl;
        Rs = R1*Rl/(R1+Rl);
        Rs0 = R1*Rs/(R1+Rs);
        if Rs0>Rl
            t = sqrt(Rs0/Rl);
        else
            t = sqrt(Rl/Rs0);
        end
        eta = ((t+1/t)/2)^2;
%         V0 = Rl/(Rs0+Rl);
        RL0   = R1; % P
        SP0   = 'P';
        Pos   = 'Head';
%         Rs0   = Rs*R1/(Rs+R1);
        V0 = Rl/(Rs0+Rl);
        f0    = ((0.8123.*Rs0./Rs).*Rs);%;
        absND = absND.*V0;
    end
    [absH] = absND;
    ND     = (funRecursionPoly(n, P));
    if isinf(eta)
        K2    = absH;
        K2(1) = 0;
    else
        Rs0   = Rs*R1/(Rs+R1);
        
        K2    = 1.27.*absH.*(Rl/(4*Rs0));%.*V0;%%;%./V0;%;%
        K2(1) = K2(1)-1;
    end
    m = length(K2);
    KK = 1i.^[0:m-1];
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
    rootsSel = rootsSel(real(rootsSel)<0);
    while length(rootsSel) > n
        rootsSel(min(abs(real(rootsSel)))==abs(real(rootsSel))) = [];
    end
    nZeroAdd = floor(nZeros/2);
    for ii=1:nZeroAdd
        rootsSel(end+1) = 0;
    end
    Fsflip = funRecursionPoly(n, rootsSel);
    NDflip = ND;
    Fs = abs(real(NDflip));
    Es = abs(real(Fsflip));
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
    km = funContinuedFractionExp(n, Z, P.*f0);
%     km(1)*km(2)
%     1.234*0.5739
%     funContinuedFractionExp(n, Z, P./1.992)
%     km = fliplr(km)./w1;
    km = (km)./w1;
%     kk = Rl/Rs;
%     km = ([1.234./kk, 0.5739.*kk]);
%     RL0 = 0.5;
    if strcmp(Pos,'Head')
        cellValueNetlist = {{'R', SP0, RL0}};
    else
        cellValueNetlist = [];
    end
    for ii=1:n
        if mod(ii, 2)
            Type = 'C';
            SP   = 'P';
        else
            Type = 'L';
            SP   = 'S';
        end
        cellValueNetlist{length(cellValueNetlist)+1} = {Type, SP, km(ii)};
    end
    if ~strcmp(Pos,'Head')
        cellValueNetlist{length(cellValueNetlist)+1} = {'R', SP0, RL0};
    end
