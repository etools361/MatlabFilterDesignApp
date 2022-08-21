%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成无向图
%--------------------------------------------------------------------------
function [h, t, e, idx, nType, n] = funCreateUndirectedGraphs(iType, Value, cellNode1, CellNode2, cellName)
[node1, node2, n] = funSimNetlistRenode(cellNode1, CellNode2);
nType = length(iType);
% 创建无向图
h   = cell(1,n);% 节点
e   = [];% 边
t   = cell(1,max(iType)+1);% 器件
idx = 0;
for ii = 1:nType
    [h, t, e, idx] = addNode(node1(ii), node2(ii), iType(ii), Value(ii), cellName(ii), h, t, e, idx);
    [h, t, e, idx] = addNode(node2(ii), node1(ii), iType(ii), Value(ii), cellName(ii), h, t, e, idx);
end

function [h, t, e, idx] = addNode(node1, node2, iType, Value, strName, h, t, e, idx)
idx = idx + 1;
e(idx).node1    = node1+1;
e(idx).node2    = node2+1;
e(idx).iType    = iType;
e(idx).Value    = Value;
e(idx).strName  = strName;
e(idx).edge     = idx;
e(idx).isMarked = 0;
h{node1+1} = [h{node1+1},idx];
t{iType+1} = [t{iType+1},idx];