%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-21(yyyy-mm-dd)
% Butterworth 滤波器综合，实现了低通原型参数计算
%--------------------------------------------------------------------------
function [cellValueNetlist, km] = funSynthesisButterworthFilter(n, Rs, Rl, fp, fs, Ap, As)
if isempty(Ap) || Ap<0
    Ap = 3;
    fprintf('Ap=%f dB\n', Ap);
end
% Butterworth参数计算
if isempty(n) || n < 2
    % 由 As计算n
    n_min = 1/log(fs/fp)*log(sqrt((10^(-0.1*As)-1)/(10^(0.1*Ap)-1)));
    n = ceil(n_min);
    fprintf('Order=%d\n', n);
else
end
epsilon = sqrt(10^(0.1*Ap)-1);
rho     = epsilon/sqrt(epsilon^2+1);
aE      = (sqrt(rho^(-2)-1))^(1/n);
if Rs==0 || Rl==0 || Rs==inf || Rl==inf
    % 一端接载
    b0 = 1;
    bm(1)   = b0*aE*2;
    for ii=2:n
        bm(ii) = 4*aE^2/bm(ii-1)*cos((ii-1)*pi/(2*n))^2;
    end
else
    % 两端接载
    if Rs>Rl
        t2 = Rs/Rl;
    else
        t2 = Rl/Rs;
    end
    kt2 = (t2-1)/(t2+1);
    aF      = (kt2*sqrt(rho^(-2)-1))^(1/n);
    b0      = aE-aF;
    b       = sqrt(4*aE*aF);
    bm      = zeros(1, n);
    bm(1)   = b0;
    for ii=2:n
        bm(ii) = 1/bm(ii-1)*(b0^2+b^2*sin((ii-1)*pi/(2*n))^2);
    end
end
am = zeros(1, n);
for ii=1:n
    am(ii) = 2*sin((2*ii-1)*pi/(2*n));
end
km         = am./bm;
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

end