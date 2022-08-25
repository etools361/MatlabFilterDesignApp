%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 滤波器综合，实现了低通原型参数计算，高通、带通、带阻转换
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisTransAndGenNetlist(fShape, TeeEn, n, Rs, Rl, fp, bw, km)
if isempty(km)
    strNetlist = [];
    return;
end
Rvs = 0;
Rvs = xor(Rvs, Rs <= Rl);
% 滤波器类型转换和网表生成
strNetlist = {};
strTemp    = sprintf('V0 V %d %d %e', n*2, 0, 1);
strNetlist = [strNetlist; strTemp];
strTemp    = sprintf('RS R %d %d %e', 1, n*2, Rs);
strNetlist = [strNetlist; strTemp];
mNode      = 1;
Rvs = xor(Rvs, ~TeeEn);
Rvs = xor(Rvs, mod(n, 2));
if Rvs
    R0 = Rl;
else
    R0 = Rs;
end
L0 = R0/(2*pi*fp);
C0 = 1/(2*pi*fp*R0);
for ii=1:n
    if Rvs
        idx = n+1-ii;
    else
        idx = ii;
    end
    odd = 1;
    if (~mod(n, 2) && (Rs<Rl || (Rs==Rl && TeeEn))) || (mod(n, 2) && TeeEn)
        odd = 0;
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
    if ~mod(ii, 2)
        odd = ~odd;% 用于产生交替信号，交替选择电容或电感
    end
    switch fShape
        case 'LPF'
            if ~odd
                strDev0 = 'L';
                Value   = km(idx)*L0;
            else
                strDev0 = 'C';
                Value   = km(idx)*C0;
            end
            mNode2  = mNode;
        case 'HPF'
            if ~odd
                strDev0 = 'C';
                Value   = 1/(km(idx))*C0;
            else
                strDev0 = 'L';
                Value   = 1/(km(idx))*L0;
            end
            mNode2  = mNode;
        case 'BPF'
            a  = bw/fp;
            if ~odd
                % L + C
                ValueL   = km(idx)/a*L0;
                ValueC   = a/km(idx)*C0;
            else
                % L // C
                ValueC   = km(idx)/a*C0;
                ValueL   = a/km(idx)*L0;
            end
            strDev0 = 'L';
            strDev  = sprintf('%s%d', strDev0, ii);
            strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), ValueL);
            strNetlist = [strNetlist; strTemp];
            strDev0 = 'C';
            mNode = mNode2;
            Value = ValueC;
        case 'BRF'
            a  = bw/fp;
            if ~odd
                % L // C
                ValueL   = km(idx)*a*L0;
                ValueC   = 1/(km(idx)*a)*C0;
            else
                % L + C
                ValueC   = km(idx)*a*C0;
                ValueL   = 1/(km(idx)*a)*L0;
            end
            if mNode(2) == 0 % 对地枝
                mNode1 = [mNode(1), n + mNode(1)];
                mNode2 = [mNode(1) + n, 0];
            else % 桥
                mNode1 = mNode;
                mNode2 = mNode;
            end
            strDev0 = 'L';
            strDev  = sprintf('%s%d', strDev0, ii);
            strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode1(1), mNode1(2), ValueL);
            strNetlist = [strNetlist; strTemp];
            strDev0 = 'C';
            Value = ValueC;
    end
    strDev  = sprintf('%s%d', strDev0, ii);
    strTemp = sprintf('%s %s %d %d %e', strDev, strDev0, mNode2(1), mNode2(2), Value);
    strNetlist = [strNetlist; strTemp];
end
netMax  = max(mNode);
strTemp = sprintf('RL R %d %d %e', netMax, 0, Rl);
strNetlist = [strNetlist; strTemp];
end