%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 重新转换节点标号
%--------------------------------------------------------------------------
function [node1, node2, n] = funSimNetlistRenode(cell_node1, cell_node2)
m = length(cell_node1);
cellAll = [cell_node1; cell_node2];
mAll    = [];
[a, b] = unique(cellAll);
n = length(a);
mAll = zeros(1, 2*m);
for ii=1:n
    ifind1 = ismember(cellAll, a(ii));
    if ~isempty(ifind1)
        mAll(ifind1) = ii-1;
    end
end
isGND = ismember(cellAll, 'GND');
if sum(isGND) == 0
    isGND   = ismember(cellAll, '0');
end
if sum(isGND) > 0
    mAll2 = mAll;
    mAll2(mAll == 0) = mAll(find(isGND==1,1));
    mAll2(isGND) = 0;
    mAll = mAll2;
else
    [a1, b1] = hist(mAll, 0:n-1);
    [a2] = find(a1 == max(a1), 1);
    b2   = b1(a2);
    if b2 ~= 0% 交换0和b2
        mAll2 = mAll;
        mAll2(mAll == 0) = b2;
        mAll2(mAll == b2) = 0;
        mAll = mAll2;
    end
end
node1 = mAll(1:m);
node2 = mAll(1+m:end);
end