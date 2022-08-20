%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 分析网表，分析网表中含有多少个电感，电流源，电压源，0欧姆电阻，和节点数
%--------------------------------------------------------------------------
function [maxNode, nL, nI, nV, nR0] = funSimNetlistAna(iType, Value, node1, node2)
nNetlist    = length(iType);
nR          = 0;
nC          = 0;
nL          = 0;
nV          = 0;
nI          = 0;
nn          = 0;
nR0         = 0;
nodeAll     = zeros(1, nNetlist);
for ii = 1:nNetlist
    % 0:V,1:I,2:R,3:L,4:C
    iTypec = iType(ii);
    if iTypec == 0
        nV = nV + 1;
    elseif iTypec == 1
        nI = nI + 1;
    elseif iTypec == 2
        if Value(ii) < 1e-18
            nR0 = nR0 + 1;
        end
        nR = nR + 1;
    elseif iTypec == 3
        nL = nL + 1;
    elseif iTypec == 4
        nC = nC + 1;
    end
    nodeAll(node1(ii)+1) = 1;
    nodeAll(node2(ii)+1) = 1;
end
nn = sum(nodeAll)-1;
maxNode = nn;
end