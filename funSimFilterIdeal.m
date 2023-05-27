%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 不同类型的滤波器仿真，理想AC特性
%--------------------------------------------------------------------------
function [IdealFreq, IdealMag, IdealPhase, P, Z, f_min] = funSimFilterIdeal(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape, f0, f1, N)
f_min      = [];
IdealFreq  = logspace(log10(f0), log10(f1), N);
s          = 1i.*2.*pi.*IdealFreq./(2.*pi.*fp);
if Rs == 0 || Rs == inf
    IdealData  = 1;
elseif Rl == 0 || Rl == inf
    IdealData  = 1;
else
    IdealData  = 1*(Rl)/(Rs+Rl);
end
PhaseTrim  = mod((n-2)*90, 360)-180;
P = [];
Z = [];
switch fShape
    case 'LPF'
        % 5:90(90), 4: 0(0), 3:-90(270), 2:-180(180)
%         s = s;
    case 'HPF'
        % 5:90(90), 4: 0(0), 3:-90(270), 2:+180(180)
        s  = 1./s;
        if PhaseTrim == -180
            PhaseTrim = 180;
        end
    case 'BPF'
        % 5:90(90), 4: 0(0), 3:-90(270), 2:+180(180)
        a  = (bw/fp);
        s = s./a+1./(s.*a);
        if PhaseTrim == -180
            PhaseTrim = 180;
        end
    case 'BRF'
        a  = (bw/fp);
        s  = a./(s+1./s);
end
switch fType % 滤波器类型
    case 'Butterworth'
        epsilon   = sqrt(10^(0.1*Ap)-1);% 通带衰减量
        aE        = (epsilon)^(-1/n);
        for ii=1:n
            k = ii;
            Zv  = aE.*exp(1i.*((2*k-1).*pi./(2*n) + pi/2));
            IdealData = IdealData.*aE./(s+Zv);
            P(ii) = Zv;
        end
        Z = inf;
    case 'Chebyshev I'
        epsilon   = sqrt(10^(0.1*Ap)-1);
        phi2      = 1/n*asinh(1/epsilon);
        IdealData = IdealData/(epsilon*2^(n-1));
        v2        = (n-1)*pi/(2*n);
        for ii=1:n
            k  = ii;
            v  = (2*k-1)*pi/(2*n);
            Zv  = (-sinh(phi2).*sin(v) + 1i.*cosh(phi2).*cos(v));
            if mod(n, 2)
                KK2 = Zv;
                C   = 1;
            else
                KK2 = sqrt(Zv.^2 + cos(v2).^2)./sin(v2);
                C   = 1./sin(v2);
            end
            P(ii) = KK2;
            IdealData = C.*IdealData.*1./(s + KK2);
%             IdealData = kk.*IdealData.*1./(w+kk.*(1.*cos(v)+1i.*1.*sin(v)));
        end
        Z = inf;
    case 'Chebyshev II' % inverse chebyshev filter
        epsilon   = 1/sqrt(10^(0.1*As)-1);% 阻带衰减量
        epsilon2  = 1/sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
        K         = cosh(1/n*acosh(epsilon2/epsilon));
%         s         = s/K;
        % 计算出阻带衰减达到最指定值时的最低频率
        v2        = (n-1)*pi/(2*n);
        if ~mod(n, 2)
            K = sqrt((K).^2 + cos(v2).^2)./sin(v2);
        end
        f_min     = fp*K;
        phi2      = 1/n*asinh(1/epsilon);
        for ii=1:n
            k  = ii;
            v  = (2*k-1)*pi/(2*n);
            Zv  = (-sinh(phi2).*sin(v) + 1i.*cosh(phi2).*cos(v));
            Zm  = -1i.*cos(v);
            if mod(n, 2)
                KK2 = K/Zv;
                KK1 = K/Zm;
                C   = 1;
            else
                KK2 = K./(sqrt((Zv).^2 + cos(v2).^2)./sin(v2));
                KK1 = sign(imag(Zm)).*K./(sqrt((Zm).^2 + cos(v2).^2)./sin(v2));
                C   = 1;
            end
            P(ii) = KK2;
            Z(ii) = KK1;
            if abs(real(KK1)) > 10*fp*2*pi
                IdealData = C.*KK2.*IdealData.*1./(s + KK2);
            else
                IdealData = C.*KK2./KK1.*IdealData.*(s + KK1)./(s + KK2);
            end
        end
    case 'Elliptic'
        % ref:Lecture Notes On Elliptic Filter Design. Sophocles J.
        es   = sqrt(10^(0.1*As)-1);% 阻带衰减量
        ep   = sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
        k1   = ep/es;
        k    = ellipdeg(n, k1);
        v2   = (n-1)/(n);
        wa   = cde(v2, k);
        wb   = 1./(k*cde((n-1)/n, k));
        aa   = wb^2;
        dd   = (-1+wb^2)/(1-wa^2);
        bb   = dd*wa^2;
        cc   = 1;
        for ii=1:n
            k0     = ii;
            u      = (2*k0-1)/n;
            KK1    = 1i./(k*cde(u, k));
            v0     = asne(1i/ep, k1)/n;
            KK2    = 1i*cde(u+v0, k);
            C = 1;
            if ~mod(n, 2)
%                 KK1 = -sign(imag(KK1)).*(sqrt(((KK1).^2 + wa.^2)./(1-wa.^2)));
%                 KK2 = (sqrt(((KK2).^2 + wa.^2)./(1-wa.^2)));
%                 KK1 = (sqrt(((KK1).^2 + wa.^2)./(1-wa.^2)));
%                 KK2 = sign(imag(KK2)).*(sqrt((wb.^2-1)./(wb.^2-(KK2).^2)));
                KK1 = sign(imag(KK1)).*sqrt((dd.*KK1.^2+bb)./(cc.*KK1.^2+aa));
                KK2 = sqrt((dd.*KK2.^2+bb)./(cc.*KK2.^2+aa));
            end
            P(ii) = KK2;
            Z(ii) = KK1;
            if abs((KK1)) > 10*2*pi*fp
                IdealData = C.*KK2.*IdealData.*1./(s + KK2);
            else
                IdealData = C.*KK2./KK1.*IdealData.*(s + KK1)./(s + KK2);
            end
        end
    case 'Bessel'
        ep   = sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
        [N2, D2, ND] = fun_bessel_thomson_polynomial(n);
        [absH] = funCalcuHjw2(ND);
        % 获取归一化频率
        [w1] = funGetBesselNormFreq(absH, ep);
        [P] = roots(fliplr(ND));
        P = P'./w1;
        Z = inf;
        IdealData = IdealData*ND(1)/w1^n;
        for ii=1:n
            IdealData = IdealData./(s - P(ii));
        end
    case 'Gaussian'
        ep   = sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
        [P, wr, absND] = fun_gaussian_polynomial(n, ep);
        P = eval(P.*wr);
        Z         = inf;
        IdealData = IdealData.*sqrt(factorial(n)).*wr.^(n);
        for ii=1:n
            IdealData = IdealData./(s - P(ii));
        end
    otherwise
        error('TBD');
end
Z(abs(Z)>100*fp*2*pi) = [];
P = [P,-real(P)+1i.*imag(P)];
Z = [Z,-real(Z)+1i.*imag(Z)];
switch fShape
    case 'LPF'
        P = P;
        Z = Z;
    case 'HPF'
        P = 1./P;
        Z = 1./Z;
    case 'BPF'
        Pa = P.*a;
        P = [(Pa+sqrt(Pa.^2-4))./2, (Pa-sqrt(Pa.^2-4))./2];
        Z = 1./Z;
    case 'BRF'
        Pa = a./P;
        P = [(Pa+sqrt(Pa.^2-4))./2, (Pa-sqrt(Pa.^2-4))./2];
        Z = [1i.*fp.*2.*pi, -1i.*fp.*2.*pi];
end
P = fp.*P;
Z = fp.*Z;
% plot(real(P), imag(P), 'x');
% hold on;
% plot(real(Z), imag(Z), 'o');grid on;xlabel('real/\delta');ylabel('imag/jw');
% hold off;
% axis square;
% xlim([-1, 1]);ylim([-1,1])
IdealMag   = 20.*log10(abs(IdealData));
IdealPhase = angle(IdealData)*180/pi + PhaseTrim;% 5:90(90), 4: 0(0), 3:-90(270), 2:-180(180)
IdealPhase = unwrap(IdealPhase, 179);
end