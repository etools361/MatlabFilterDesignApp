%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 不同类型的滤波器综合
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape)
switch fType % 滤波器类型
    case 'Butterworth'
        [cellValueNetlist, km] = funSynthesisButterworthFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Chebyshev I'
        [cellValueNetlist, km] = funSynthesisChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Chebyshev II'
        [cellValueNetlist, km] = funSynthesisInverseChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Elliptic'
        [cellValueNetlist, km] = funSynthesisEllipticFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Bessel'
        [cellValueNetlist, km] = funSynthesisBesselFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Gaussian'
        [cellValueNetlist, km] = funSynthesisGaussianFilter(n, Rs, Rl, fp, fs, Ap, As);
    otherwise
        error('TBD');
        km = [];
end
% [strNetlist] = funSynthesisTransAndGenNetlist(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, km);
[strNetlist] = funSynthesisTransAndGenNetlist2(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, cellValueNetlist);
end