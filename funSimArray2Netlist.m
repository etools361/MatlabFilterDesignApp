%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 将网表转换为数组形式
%--------------------------------------------------------------------------
function [strNetlist] = funSimArray2Netlist(iType, Value, node1, node2, cellName)
n    = length(iType);
strNetlist = [];
for ii = 1:n
    [Type] = funSimiType2Type(iType(ii));
    strNetlist{ii} = sprintf('%s %s %s %s %e\n' , cellName{ii}, Type, node1{ii}, node2{ii}, Value(ii));
end
end