%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 滤波器综合，实现了低通原型参数计算，高通、带通、带阻转换
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisTransAndGenNetlist(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, km)
if isempty(km)
    strNetlist = [];
    return;
end
km2  = 0;
km2c = 1;
switch fType % 滤波器类型
    case 'Butterworth'
        isZeros = 0;
    case 'Chebyshev I'
        isZeros = 0;
    case 'Chebyshev II'
        isZeros = 1;
        km2 = km(3:3:end);
        km(3:3:end)  = [];
    otherwise
        error('TBD');
        km = [];
end
mkm2c = length(km2);
Rvs = 0;
% 滤波器类型转换和网表生成
strNetlist = {};
if Rs == inf
    strTemp    = sprintf('I0 I %d %d %e', n*20+1, 0, 1/Rl);
    strNetlist = [strNetlist; strTemp];
    strTemp    = sprintf('RS R %d %d %e', 1, n*20+1, 0);
    strNetlist = [strNetlist; strTemp];
else
    if Rl == 0
        strTemp    = sprintf('V0 V %d %d %e', n*20+1, 0, Rs);
    else
        strTemp    = sprintf('V0 V %d %d %e', n*20+1, 0, 1);
    end
    strNetlist = [strNetlist; strTemp];
    strTemp    = sprintf('RS R %d %d %e', 1, n*20+1, Rs);
    strNetlist = [strNetlist; strTemp];
end
mNode      = 1;
OneEndPort = 0;
if Rl==0 || Rl==inf
    Rvs = 1;
    OneEndPort = 1;
    Rvs = xor(Rvs, ~mod(n, 2));
elseif Rs==0 || Rs==inf
    Rvs = 0;
    OneEndPort = 1;
    Rvs = xor(Rvs, ~mod(n, 2));
else
    Rvs = xor(Rvs, Rs <= Rl);
    Rvs = xor(Rvs, ~TeeEn);
end
Rvs = xor(Rvs, mod(n, 2));
if Rvs
    R0 = Rl;
else
    R0 = Rs;
end
if Rl==0 || Rl==inf
    R0 = Rs;
elseif Rs==0 || Rs==inf
    R0 = Rl;
end
L0 = R0/(2*pi*fp);
C0 = 1/(2*pi*fp*R0);
for ii=1:n
    if Rvs
        idx   = n+1-ii;
        km2cx = mkm2c+1-km2c;
    else
        idx = ii;
        km2cx = km2c;
    end
    odd = 1;
    if OneEndPort
        if Rl==0 || Rs==0 % T
            odd = 0;
%             odd = xor(odd, ~mod(n, 2));
        elseif Rl==inf || Rs==inf % PI
            odd = 1;
        end
        if Rs == 0 || Rs == inf
        else
            odd = xor(odd, ~mod(n, 2));
        end
    else
        if (~mod(n, 2) && (Rs<Rl || (Rs==Rl && TeeEn))) || (mod(n, 2) && TeeEn)
            odd = 0;
        end
    end
    if mod(ii+odd, 2)
        mNode   = [ceil(ii/2), ceil(ii/2)+1];
        if ~odd
            mNode1   = [ii, ii+1];
            mNode2   = [ii+1, ii+2];
        else
            mNode1   = [ii-1, ii];
            mNode2   = [ii, ii+1];
        end
    else
        mNode   = [max(mNode), 0];
        mNode1  = mNode;
        mNode2  = mNode;
    end
    odd = xor(odd, ~mod(ii, 2));
    switch fShape
        case 'LPF'
            if ~odd
                strDev0 = 'L';
                Value   = km(idx)*L0;Value2 = km2(km2cx)*C0;
            else
                strDev0 = 'C';
                Value   = km(idx)*C0;Value2 = km2(km2cx)*L0;
            end
            mNode2  = mNode;
        case 'HPF'
            if ~odd
                strDev0 = 'C';
                Value   = 1/(km(idx))*C0;Value2 = 1/(km2(km2cx))*L0;
            else
                strDev0 = 'L';
                Value   = 1/(km(idx))*L0;Value2 = 1/(km2(km2cx))*C0;
            end
            mNode2  = mNode;
        case 'BPF'
            a  = bw/fp;
            if isZeros && ii ~= n && ~mod(ii, 2)
                % 带0点的臂综合，需要频率转换
                W          = 1/((km(idx)/a)*(km2(km2cx)/a));
                Beta1      = 1+W/2+sqrt(W+1/4*W^2);
                L2         = 1/((km2(km2cx)/a)*(1+Beta1));
                C2         = 1/(Beta1*L2);
                km(idx)    = L2;
                km2(km2cx) = C2;
                a = 1;
            end
            if ~odd
                % L + C
                ValueL   = km(idx)/a*L0;Value2L = km2(km2cx)/a*C0;
                ValueC   = a/km(idx)*C0;Value2C = a/km2(km2cx)*L0;
            else
                % L // C
                ValueC   = km(idx)/a*C0;Value2C = km2(km2cx)/a*L0;
                ValueL   = a/km(idx)*L0;Value2L = a/km2(km2cx)*C0;
            end
            strDev0 = 'L';
            strDev  = sprintf('%s%d', strDev0, ii);
            if TeeEn && isZeros && ii ~= n && ~mod(ii, 2)
                mNode2_Mid1 = n*3+mNode1(1);
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode2_Mid1, ValueL);
            else
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), ValueL);
            end
            strNetlist = [strNetlist; strTemp];
            if isZeros && ii ~= n && ~mod(ii, 2)
                if strDev0 == 'L'
                    strDev0 = 'C';
                else
                    strDev0 = 'L';
                end
                strDev  = sprintf('%s%d', strDev0, ii);
                if TeeEn
                    strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2_Mid1, mNode1(2), Value2L);
                else
                    strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), Value2L);
                end
                    strNetlist = [strNetlist; strTemp];
            end
            strDev0 = 'C';
            mNode = mNode2;
            Value = ValueC;Value2 = Value2C;
        case 'BRF'
            a  = bw/fp;
            if isZeros && ii ~= n && ~mod(ii, 2)
                % 带0点的臂综合，需要频率转换
                W          = 1/(1/(km(idx)*a)*1/(a*km2(km2cx)));
                Beta1      = 1+W/2+sqrt(W+1/4*W^2);
                L2         = 1/(1/(a*km(idx))*(1+Beta1));
                C2         = 1/(Beta1*L2);
                km(idx)    = L2;
                km2(km2cx) = C2;
                a = 1;
            end
            if ~odd
                % L // C
                ValueL   = km(idx)*a*L0;     Value2L = km2(km2cx)*a*C0;
                ValueC   = 1/(km(idx)*a)*C0; Value2C = 1/(km2(km2cx)*a)*L0;
            else
                % L + C
                ValueC   = km(idx)*a*C0;     Value2C = km2(km2cx)*a*L0;
                ValueL   = 1/(km(idx)*a)*L0; Value2L = 1/(km2(km2cx)*a)*C0;
            end
            if mNode(2) == 0 % 对地枝
                if TeeEn && isZeros && ii ~= n && ~mod(ii, 2)
                    mNode1 = [mNode(1), 0];
                    mNode2 = [mNode(1), 0];
                else
                    mNode1 = [mNode(1), 2*n + mNode(1)];
                    mNode2 = [mNode(1) + 2*n, 0];
                end
            else % 桥
                mNode1 = mNode;
                if isZeros && ii ~= n && ~mod(ii, 2)
                else
                    mNode2 = mNode;
                end
            end
            strDev0 = 'L';
            strDev  = sprintf('%s%d', strDev0, ii);
            if TeeEn && isZeros && ii ~= n && ~mod(ii, 2)
                mNode2_Mid1 = n*3+mNode1(1);
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode2_Mid1, ValueL);
            else
                strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), ValueL);
            end
            strNetlist = [strNetlist; strTemp];
            if isZeros && ii ~= n && ~mod(ii, 2)
                if strDev0 == 'L'
                    strDev0 = 'C';
                else
                    strDev0 = 'L';
                end
                strDev  = sprintf('%s%d', strDev0, ii);
                if TeeEn
                    strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2_Mid1, mNode1(2), Value2L);
                else
                    strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), Value2L);
                end
                strNetlist = [strNetlist; strTemp];
                mNode = mNode2;
            end
            strDev0 = 'C';
            Value = ValueC; Value2 = Value2C;
    end
    strDev  = sprintf('%s%d', strDev0, ii);
    if TeeEn && isZeros && ii ~= n && ~mod(ii, 2)
        mNode2_Mid2 = n*4+mNode1(1);
        strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2(1), mNode2_Mid2, Value);
    else
        strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2(1), mNode2(2), Value);
    end
    strNetlist = [strNetlist; strTemp];
    if isZeros && ii ~= n && ~mod(ii, 2)
        if strDev0 == 'L'
            strDev0 = 'C';
        else
            strDev0 = 'L';
        end
        strDev  = sprintf('%s%d', strDev0, ii);
        if TeeEn
            strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2_Mid2, mNode2(2), Value2);
        else
            strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2(1), mNode2(2), Value2);
        end
        km2c    = km2c + 1;
        if km2c > mkm2c
            km2c = mkm2c;
        end
        strNetlist = [strNetlist; strTemp];
    end
end
netMax  = max(mNode);
strTemp = sprintf('RL R %d %d %e', netMax, 0, Rl);
strNetlist = [strNetlist; strTemp];
end