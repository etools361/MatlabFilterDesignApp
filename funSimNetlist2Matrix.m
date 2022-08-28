%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 由网表转换为矩阵
%--------------------------------------------------------------------------
function [M, N, MV, MX] = funSimNetlist2Matrix(iType, Value, node1, node2, maxNode, nL, nI, nV, nR0, cellName)
nNetlist    = length(iType);
% | MC 0  |   d       | MG  MD |
% |       | * -  MX + |        | * MX = MV
% | 0  ML |   dt      | MD' MI |
% 
% nR0表示R=0的数量
% nL表示电感数量
% nV表示电压源数量
% nL表示电流源数量
% maxNode表示节点数量-1
MC     = zeros(maxNode, maxNode);
MG     = zeros(maxNode, maxNode);
nVI    = nV+nI+nR0+nL;
% NVI    = zeros(nVI, nVI);
ML     = zeros(nVI, nVI);
MI     = zeros(nVI, nVI);
MD     = zeros(maxNode, nVI);
MX     = cell(1, nVI+maxNode);
MV     = zeros(1, nVI+maxNode);
iMD    = 0;
for ii=1:maxNode
    MX{ii} = sprintf('V%dn', ii);
end
for ii = 1:nNetlist
    % 0:V,1:I,2:R,3:L,4:C
    iTypec = iType(ii);
    if iTypec == 0 % V
        iMD = iMD + 1;
        [ML, MD] = funNode2MatR(node1(ii), node2(ii), ML, MD, iMD, 0);%node1, node2, MI, MD, iMD, Value
        MV(iMD + maxNode) = Value(ii);
        MX{iMD + maxNode} = sprintf('i%s', cellName{ii});
    elseif iTypec == 1 % I
        iMD = iMD + 1;
        [MI, MD] = funNode2MatR(node1(ii), node2(ii), MI, MD, iMD, 1);
%         [MI, MD] = funNode2MatR(node1(ii), node2(ii), MI, MD, iMD, Value(ii));
        MV(iMD + maxNode) = Value(ii);
        MX{iMD + maxNode} = sprintf('i%s', cellName{ii});
    elseif iTypec == 2 % R
        if Value(ii)>1e-18
            [MG] = funNode2Mat(node1(ii), node2(ii), MG, 1/Value(ii));
        else
            iMD = iMD + 1;
            [MI, MD] = funNode2MatR(node1(ii), node2(ii), MI, MD, iMD, 0);%node1, node2, MI, MD, iMD, Value
            MX{iMD + maxNode} = sprintf('i%s', cellName{ii});
        end
    elseif iTypec == 3 % L
        iMD = iMD + 1;
        [ML, MD] = funNode2MatR(node1(ii), node2(ii), ML, MD, iMD, Value(ii));
        MX{iMD + maxNode} = sprintf('i%s', cellName{ii});
    elseif iTypec == 4 % C
        [MC] = funNode2Mat(node1(ii), node2(ii), MC, Value(ii));
    end
end
% MC, MG, ML, MI, MD
% | MC 0  |   d       | MG  MD |
% |       | * -  MX + |        | * MX = MV
% | 0  ML |   dt      | MD' MI |
M = [MC, zeros(maxNode, nVI);zeros(maxNode, nVI)', ML];
MDT = -MD';
% 补丁，电流源矩阵修正
[a, ~] = find(MI==1);
MDT(a,:) = 0;
N = [MG, MD; MDT, MI];
end

function [MI] = funNode2Mat(node1, node2, MI, Value)
    if node1 == 0
        MI(node2, node2) = MI(node2, node2) + Value;
    elseif node2 == 0
        MI(node1, node1) = MI(node1, node1) + Value;
    else
        MI(node1, node1) = MI(node1, node1) + Value;
        MI(node2, node2) = MI(node2, node2) + Value;
        MI(node1, node2) = MI(node1, node2) - Value;
        MI(node2, node1) = MI(node2, node1) - Value;
    end
end

function [MI, MD] = funNode2MatR(node1, node2, MI, MD, iMD, Value)
    MI(iMD, iMD) = Value;
    if node1 == 0
        MD(node2, iMD) = 1;
    elseif node2 == 0
        MD(node1, iMD) = -1;
    else
        MD(node1, iMD) = -1;
        MD(node2, iMD) = 1;
    end
end