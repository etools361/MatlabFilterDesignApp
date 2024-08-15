%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 不同类型的滤波器综合
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, Apr, Asr, bw, fShape)
switch fType % 滤波器类型
    case 'Butterworth'
        [cellValueNetlist, km] = funSynthesisButterworthFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    case 'Chebyshev I'
        [cellValueNetlist, km] = funSynthesisChebyshevFilter(n, Rs, Rl, fp, fs, Apr, Asr);
        v2  = (n-1)*pi/(2*n);
        fpx = cos(1/n*acos(sqrt(10^(Ap/10)-1)/sqrt(10^(Apr/10)-1)));
        if ~mod(n,2)
            fpx = sqrt(fpx^2-cos(v2)^2)/sin(v2);
        end
        fp  = fp/fpx;
    case 'Chebyshev II'
        [cellValueNetlist, km] = funSynthesisInverseChebyshevFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    case 'Elliptic'
        [cellValueNetlist, km] = funSynthesisEllipticFilter(n, Rs, Rl, fp, fs, Apr, Asr);
        es   = sqrt(10^(0.1*Asr)-1);% 阻带衰减量
        ep   = sqrt(10^(0.1*Ap)-1);% 截止频率处衰减量
        epr  = sqrt(10^(0.1*Apr)-1);% 截止频率处衰减量
        k1   = epr/es;
        k    = ellipdeg(n, k1);
        fpx = cde(1/n*acde(ep/epr,k1),k);
        if ~mod(n,2)
            v2   = (n-1)/(n);
            wa   = cde(v2, k);
            wb   = 1./(k*wa);
            aa   = wb^2;
            dd   = (-1+wb^2)/(1-wa^2);
            bb   = dd*wa^2;
            cc   = 1;
            fpx = sqrt((-bb+dd*(fpx)^2)/(-cc*(fpx)^2+aa));
        end
        fp  = fp/fpx;
    case 'Bessel'
        [cellValueNetlist, km] = funSynthesisBesselFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    case 'Gaussian'
        [cellValueNetlist, km] = funSynthesisGaussianFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    case 'Legendre'
        [cellValueNetlist, km] = funSynthesisLegendreFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    case 'LinearAmp'
        addpath('../chebfun-master');
        [cellValueNetlist, km, Rs, Rl] = funSynthesisLinearAmpFilter(n, Rs, Rl, fp, fs, Ap, Asr);
    otherwise
        error('TBD');
        km = [];
end
% [strNetlist] = funSynthesisTransAndGenNetlist(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, km);
[strNetlist] = funSynthesisTransAndGenNetlist2(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, cellValueNetlist);
end