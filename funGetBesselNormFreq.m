%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-02-12(yyyy-mm-dd)
% 计算贝塞尔滤波器的归一化频率
%--------------------------------------------------------------------------
function [w1] = funGetBesselNormFreq(absH, epsilon)
r3dB    = absH;
r3dB(1) = r3dB(1)*(-epsilon^2);
[rr3dB] = roots(fliplr(r3dB));
imag0   = rr3dB(abs(imag(rr3dB))<1e-6);
w1      = imag0(imag0>0);
end