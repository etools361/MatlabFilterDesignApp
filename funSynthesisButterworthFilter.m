%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-21(yyyy-mm-dd)
% Butterworth 滤波器综合
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisButterworthFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As)
if isempty(Ap) || Ap>=0
    Ap = -3;
    fprintf('Ap=%f dB\n', Ap);
end
if isempty(n) || n < 2
    % 由 As计算n
    n_min = 1/log(fs/fp)*log(sqrt((10^(-0.1*As)-1)/(10^(-0.1*Ap)-1)));
    n = ceil(n_min);
    fprintf('Order=%d\n', n);
else
end
if Rs==0
elseif Rl==0
elseif Rs==inf
elseif Rl==inf
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
aE = (sqrt(rho^(-2)-1))^(1/n);
b0 = aE-aF;
b  = sqrt(4*aE*aF);
bm = zeros(1, n);
bm(1) = b0;
for ii=2:n
    bm(ii) = 1/bm(ii-1)*(b0^2+b^2*sin((ii-1)*pi/(2*n))^2);
end
am = zeros(1, n);
for ii=1:n
    am(ii) = 2*sin((2*ii-1)*pi/(2*n));
end
km = am./bm;
strNetlist = {};
strTemp = sprintf('V0 V %d %d %e', n+1, 0, 1);
strNetlist = [strNetlist; strTemp];
strTemp = sprintf('RS R %d %d %e', 1, n+1, Rs);
strNetlist = [strNetlist; strTemp];
netMax = 0;
R0 = Rs;
L0 = R0/(2*pi*fp);
C0 = 1/(2*pi*fp*R0);
switch fType
    case 'LPF'
        if TeeEn
            % Tee结构
            for ii=1:n
                if mod(ii, 2) == 1
                    strDev0 = 'L';
                    mNode   = [ceil(ii/2), ceil(ii/2)+1];
                    if mod(n, 2)
                        Value   = km(ii)*L0;
                    else
                        Value   = km(n+1-ii)*(Rl/Rs)*L0;
                    end
                else
                    strDev0 = 'C';
                    mNode   = [ceil(ii/2)+1, 0];
                    if mod(n, 2)
                        Value   = km(ii)*C0;
                    else
                        Value   = km(n+1-ii)/(Rl/Rs)*C0;
                    end
                end
                strDev  = sprintf('%s%d', strDev0, ii);
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode(1), mNode(2), Value);
                strNetlist = [strNetlist; strTemp];
            end
        else
            % PI结构
            for ii=1:n
                if mod(ii, 2) == 1
                    strDev0 = 'C';
                    mNode   = [ceil(ii/2), 0];
                    if mod(n, 2)
                        Value   = km(n+1-ii)/(Rl/Rs)*C0;
                    else
                        Value   = km(n+1-ii)*(Rl/Rs)*C0;
                    end
                else
                    strDev0 = 'L';
                    mNode   = [ceil(ii/2), ceil(ii/2)+1];
                    if mod(n, 2)
                        Value   = km(n+1-ii)*(Rl/Rs)*L0;
                    else
                        Value   = km(n+1-ii)/(Rl/Rs)*L0;
                    end
                end
                strDev  = sprintf('%s%d', strDev0, ii);
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode(1), mNode(2), Value);
                strNetlist = [strNetlist; strTemp];
            end
        end
    case 'HPF'
    case 'BPF'
    case 'BRF'
end
netMax  = max(mNode);
strTemp = sprintf('RL R %d %d %e', netMax, 0, Rl);
strNetlist = [strNetlist; strTemp];
end