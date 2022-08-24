%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-21(yyyy-mm-dd)
% Butterworth 滤波器综合，实现了低通原型参数计算
%--------------------------------------------------------------------------
function [km] = funSynthesisButterworthFilter(n, Rs, Rl, fp, fs, Ap, As)
if isempty(Ap) || Ap>=0
    Ap = -3;
    fprintf('Ap=%f dB\n', Ap);
end
% Butterworth参数计算
if isempty(n) || n < 2
    % 由 As计算n
    n_min = 1/log(fs/fp)*log(sqrt((10^(-0.1*As)-1)/(10^(-0.1*Ap)-1)));
    n = ceil(n_min);
    fprintf('Order=%d\n', n);
else
end
if Rs==0
    % TBD
    warning('TBD');
elseif Rl==0
    % TBD
    warning('TBD');
elseif Rs==inf
    % TBD
    warning('TBD');
elseif Rl==inf
    % TBD
    warning('TBD');
else
    if Rs>Rl
        t2 = Rs/Rl;
    else
        t2 = Rl/Rs;
    end
end
epsilon = sqrt(10^(-0.1*Ap)-1);
rho     = epsilon/sqrt(epsilon^2+1);
aF      = ((t2-1)/(t2+1)*sqrt(rho^(-2)-1))^(1/n);
aE      = (sqrt(rho^(-2)-1))^(1/n);
b0      = aE-aF;
b       = sqrt(4*aE*aF);
bm      = zeros(1, n);
bm(1)   = b0;
for ii=2:n
    bm(ii) = 1/bm(ii-1)*(b0^2+b^2*sin((ii-1)*pi/(2*n))^2);
end
am = zeros(1, n);
for ii=1:n
    am(ii) = 2*sin((2*ii-1)*pi/(2*n));
end
km         = am./bm;

end