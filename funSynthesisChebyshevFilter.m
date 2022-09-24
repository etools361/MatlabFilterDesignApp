%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-09-02(yyyy-mm-dd)
% Chebyshev 滤波器综合，实现了低通原型参数的计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km] = funSynthesisChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As)
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
    % calcu Fs
    epsilon   = sqrt(10^(0.1*Ap)-1);
    phi2      = 1/n*asinh(1/epsilon);
    B         = 1-((t^2-1)/(t^2+1))^2;
    m         = (1-B)/epsilon^2;
    h         = (sqrt(1+m)+sqrt(m))^(1/n);
    Zv        = zeros(1, n);
    Rv        = zeros(1, n);
    for ii=1:n
        k  = ii;
        v  = (2*k-1)*pi/(2*n);
        Zv(ii)  = (-sinh(phi2).*sin(v) + 1i.*cosh(phi2).*cos(v)); % 8.52
        Rv(ii)  = 0.5*(h*exp(1i*(pi/2+v))+1/h*exp(1i*(pi/2-v)));
    end
    if ~mod(n, 2)
        v2 = (n-1)*pi/(2*n);
        Zv1 = -sqrt(Zv.^2+cos(v2).^2)./sin(v2);
        Rv1 = -sqrt(Rv.^2+cos(v2).^2)./sin(v2);
    else
        Zv1 = Zv;
        Rv1 = Rv;
    end
    Fs = funRecursionPoly(n, Zv1);
    Es = funRecursionPoly(n, Rv1);
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


