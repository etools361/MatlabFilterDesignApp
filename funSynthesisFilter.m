%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 不同类型的滤波器综合
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape)
switch fType % 滤波器类型
    case 'Butterworth'
        [km] = funSynthesisButterworthFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Chebyshev I'
        [km] = funSynthesisChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As);
    case 'Chebyshev II'
        [km] = funSynthesisInverseChebyshevFilter(n, Rs, Rl, fp, fs, Ap, As);
    otherwise
        error('TBD');
        km = [];
end
[strNetlist] = funSynthesisTransAndGenNetlist(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, km);
end