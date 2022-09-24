%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-25(yyyy-mm-dd)
% 滤波器综合，实现了低通原型参数计算，高通、带通、带阻转换
%--------------------------------------------------------------------------
function [strNetlist] = funSynthesisTransAndGenNetlist2(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, cellValueNetlist)
if Rs == inf
    HeadNetlist = {{'I', 'P', 1/Rl},{'R', 'S', 0}};
else
    if Rl == 0
        HeadNetlist = {{'V', 'P', Rs},{'R', 'S', Rs}};
    else
        HeadNetlist = {{'V', 'P', 1},{'R', 'S', Rs}};
    end
end
TailNetlist = {{'R', 'P', Rl}};
nDev = length(cellValueNetlist);
% --------------------Rs and Rl---------------------------
Rvs = 0;
Rvs = xor(Rvs, ~TeeEn);
Rvs = xor(Rvs, Rs <= Rl);
Rvs = xor(Rvs, mod(n, 2));
if Rl==0
    R0    = Rs;
    TeeEn = 1;
    cellValueNetlist = fliplr(cellValueNetlist);
elseif Rl==inf
    R0    = Rs;
    TeeEn = 0;
    cellValueNetlist = fliplr(cellValueNetlist);
elseif Rs==0
    R0    = Rl;
    TeeEn = 1;
elseif Rs==inf
    R0    = Rl;
    TeeEn = 0;
else
    if Rvs
        R0 = Rl;
    else
        R0 = Rs;
        cellValueNetlist = fliplr(cellValueNetlist);
    end
end
% --------------------T or PI---------------------------
if TeeEn
    for ii=1:nDev
        Dev = cellValueNetlist{ii};
        Dev = funAltDev(Dev); 
        cellValueNetlist{ii} = Dev;
    end
end

L0 = R0/(2*pi*fp);
C0 = 1/(2*pi*fp*R0);
% --------------------Frequency Response---------------------------
cellValueNetlistMain = cellValueNetlist;
switch fShape
    case 'LPF'
        for ii=1:nDev
            Dev = cellValueNetlistMain{ii};
            switch Dev{1}
                case 'L'
                    cellValueNetlistMain{ii}{3} = cellValueNetlistMain{ii}{3}*L0;
                case 'C'
                    cellValueNetlistMain{ii}{3} = cellValueNetlistMain{ii}{3}*C0;
                case {'LCS', 'LCP'}
                    cellValueNetlistMain{ii}{3} = cellValueNetlistMain{ii}{3}*L0;
                    cellValueNetlistMain{ii}{4} = cellValueNetlistMain{ii}{4}*C0;
            end
        end
    case 'HPF'
        for ii=1:nDev
            Dev = cellValueNetlistMain{ii};
            switch Dev{1}
                case 'L'
                    cellValueNetlistMain{ii}{3} = 1/cellValueNetlistMain{ii}{3}*C0;
                    cellValueNetlistMain{ii}{1} = 'C';
                case 'C'
                    cellValueNetlistMain{ii}{3} = 1/cellValueNetlistMain{ii}{3}*L0;
                    cellValueNetlistMain{ii}{1} = 'L';
                case {'LCS', 'LCP'}
                    LSP = cellValueNetlistMain{ii}{3};
                    CSP = cellValueNetlistMain{ii}{4};
                    cellValueNetlistMain{ii}{3} = 1/CSP*L0;
                    cellValueNetlistMain{ii}{4} = 1/LSP*C0;
            end
        end
    case 'BPF'
        for ii=1:nDev
            Dev = cellValueNetlistMain{ii};
            a  = bw/fp;
            Temp = cellValueNetlistMain{ii}{3};
            switch Dev{1}
                case 'L'
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = Temp/a*L0;
                            cellValueNetlistMain{ii}{4} = a/Temp*C0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = a/Temp*C0;
                            cellValueNetlistMain{ii}{4} = Temp/a*L0;
                    end
                    cellValueNetlistMain{ii}{1} = 'LCS';
                case 'C'
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = a/Temp*C0;
                            cellValueNetlistMain{ii}{4} = Temp/a*L0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = a/Temp*L0;
                            cellValueNetlistMain{ii}{4} = Temp/a*C0;
                    end
                    cellValueNetlistMain{ii}{1} = 'LCP';
                case {'LCS', 'LCP'}
%                     % 带0点的臂综合，需要频率转换
                    switch Dev{2}
                        case 'S'
                            LSP = cellValueNetlistMain{ii}{3};
                            CSP = cellValueNetlistMain{ii}{4};
                        case 'P'
                            CSP = cellValueNetlistMain{ii}{3};
                            LSP = cellValueNetlistMain{ii}{4};
                    end
                    W          = 1/((LSP/a)*(CSP/a));
                    Beta1      = 1+W/2+sqrt(W+1/4*W^2);
                    L2         = 1/((CSP/a)*(1+Beta1));
                    C2         = 1/(Beta1*L2);
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = L2*L0;
                            cellValueNetlistMain{ii}{4} = C2*C0;
                            cellValueNetlistMain{ii}{5} = 1/C2*L0;
                            cellValueNetlistMain{ii}{6} = 1/L2*C0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = 1/L2*L0;
                            cellValueNetlistMain{ii}{4} = 1/C2*C0;
                            cellValueNetlistMain{ii}{5} = C2*L0;
                            cellValueNetlistMain{ii}{6} = L2*C0;
                    end
                    cellValueNetlistMain{ii}{1} = sprintf('%s%s%s', Dev{1}(1:2), Dev{1}(1:2), Dev{1}(3));
            end
        end
    case 'BRF'
        for ii=1:nDev
            Dev = cellValueNetlistMain{ii};
            a  = bw/fp;
            Temp = cellValueNetlistMain{ii}{3};
            switch Dev{1}
                case 'L'
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = Temp/a*L0;
                            cellValueNetlistMain{ii}{4} = a/Temp*C0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = a/Temp*C0;
                            cellValueNetlistMain{ii}{4} = Temp/a*L0;
                    end
                    cellValueNetlistMain{ii}{1} = 'LCP';
                case 'C'
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = a/Temp*C0;
                            cellValueNetlistMain{ii}{4} = Temp/a*L0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = a/Temp*L0;
                            cellValueNetlistMain{ii}{4} = Temp/a*C0;
                    end
                    cellValueNetlistMain{ii}{1} = 'LCS';
                case {'LCS', 'LCP'}
%                     % 带0点的臂综合，需要频率转换
                    switch Dev{2}
                        case 'S'
                            LSP = cellValueNetlistMain{ii}{3};
                            CSP = cellValueNetlistMain{ii}{4};
                        case 'P'
                            CSP = cellValueNetlistMain{ii}{3};
                            LSP = cellValueNetlistMain{ii}{4};
                    end
                    W          = 1/(1/(LSP/a)*1/(CSP/a));
                    Beta1      = 1+W/2+sqrt(W+1/4*W^2);
                    L2         = 1/(1/(LSP/a)*(1+Beta1));
                    C2         = 1/(Beta1*L2);
                    switch Dev{2}
                        case 'S'
                            cellValueNetlistMain{ii}{3} = L2*L0;
                            cellValueNetlistMain{ii}{4} = C2*C0;
                            cellValueNetlistMain{ii}{5} = 1/C2*L0;
                            cellValueNetlistMain{ii}{6} = 1/L2*C0;
                        case 'P'
                            cellValueNetlistMain{ii}{3} = 1/L2*L0;
                            cellValueNetlistMain{ii}{4} = 1/C2*C0;
                            cellValueNetlistMain{ii}{5} = C2*L0;
                            cellValueNetlistMain{ii}{6} = L2*C0;
                    end
                    cellValueNetlistMain{ii}{1} = sprintf('%s%s%s', Dev{1}(1:2), Dev{1}(1:2), Dev{1}(3));
            end
        end
end
cellValueNetlist = [HeadNetlist, cellValueNetlistMain, TailNetlist];
[strNetlist] = funLPF2NormNetlist(cellValueNetlist);

function [Dev] = funAltDev(Dev)
M = {'L', 'C', 'LCP', 'LCS'};
N = {'C', 'L', 'LCS', 'LCP'};
[a, b] = ismember(Dev{1}, M);
if a
    Dev{1} = N{b};
    if b == 3 || b == 4
        Dev(3:end) = fliplr(Dev(3:end));
    end
end
M = {'S', 'P'};
N = {'P', 'S'};
[a, b] = ismember(Dev{2}, M);
if a
    Dev{2} = N{b};
end


function [strNetlist] = funLPF2NormNetlist(cellValueNetlist)
strNetlist = [];
m = length(cellValueNetlist);
MainNodeMax = m;
Node    = [1, 0];
NodeCnt = 1;
for ii=1:m
    cellValue = cellValueNetlist{ii};
    Type   = cellValue{1};
    SP     = cellValue{2};
    mValue = cellValue(3:end);
    switch SP
        case 'S'
            Node = [NodeCnt, NodeCnt+1];
            NodeCnt = NodeCnt + 1;
        case 'P'
            Node = [NodeCnt, 0];
            NodeCnt = NodeCnt;
    end
    switch Type
        case {'V', 'I'}
            if ii == 1
                strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type, 0, Type, Node(1), Node(2), mValue{1})}];
            else
                strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type, ii, Type, Node(1), Node(2), mValue{1})}];
            end
        case 'R'
            if ii == 2
                strNetlist = [strNetlist; {sprintf('RS %s %d %d %e', Type, Node(1), Node(2), mValue{1})}];
            elseif ii == m
                strNetlist = [strNetlist; {sprintf('RL %s %d %d %e', Type, Node(1), Node(2), mValue{1})}];
            else
                strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type, ii, Type, Node(1), Node(2), mValue{1})}];
            end
        case {'L', 'C'}
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type, ii, Type, Node(1), Node(2), mValue{1})}];
        case {'RLP', 'RCP', 'LCP'}
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),             ii, Type(1), Node(1), Node(2), mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+MainNodeMax, Type(2), Node(1), Node(2), mValue{2})}];
        case {'RLS', 'RCS', 'LCS'}
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),             ii, Type(1), Node(1), ii+MainNodeMax, mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+MainNodeMax, Type(2), ii+MainNodeMax, Node(2), mValue{2})}];
        case 'RLCP'
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),               ii, Type(1), Node(1), Node(2), mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+1*MainNodeMax, Type(2), Node(1), Node(2), mValue{2})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(3), ii+2*MainNodeMax, Type(3), Node(1), Node(2), mValue{3})}];
        case 'RLCS'
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),               ii, Type(1),               Node(1), ii+1*MainNodeMax, mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+1*MainNodeMax, Type(2), ii+1*MainNodeMax, ii+2*MainNodeMax, mValue{2})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(3), ii+2*MainNodeMax, Type(3), ii+2*MainNodeMax, Node(2),               mValue{3})}];
        case 'LCLCP'
            %   o--+---L1---+---L3----+---o
            %      |        |         |
            %      +---C2---+---C4----+
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),               ii, Type(1),               Node(1), ii+1*MainNodeMax, mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+1*MainNodeMax, Type(2),               Node(1), ii+1*MainNodeMax, mValue{2})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(3), ii+2*MainNodeMax, Type(3), ii+1*MainNodeMax, Node(2),               mValue{3})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(4), ii+3*MainNodeMax, Type(4), ii+1*MainNodeMax, Node(2),               mValue{4})}];
        case 'LCLCS'
            %   o--+---L1------C2----+---o
            %      |                 |
            %      +---L3------C4----+
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(1),               ii, Type(1), Node(1),               ii+1*MainNodeMax, mValue{1})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(2), ii+1*MainNodeMax, Type(2), ii+1*MainNodeMax, Node(2)              , mValue{2})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(3), ii+2*MainNodeMax, Type(3), Node(1),               ii+2*MainNodeMax, mValue{3})}];
            strNetlist = [strNetlist; {sprintf('%s%d %s %d %d %e', Type(4), ii+3*MainNodeMax, Type(4), ii+2*MainNodeMax, Node(2),               mValue{4})}];
    end
end

