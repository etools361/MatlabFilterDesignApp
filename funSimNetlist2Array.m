%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 将网表转换为数组形式
%--------------------------------------------------------------------------
function [iType, Value, node1, node2, cellName] = funSimNetlist2Array(strNetlist)
nNetlist    = length(strNetlist);
cellNetlist = cell(nNetlist, 1);
cellName    = cell(nNetlist, 1);
iType       = zeros(nNetlist, 1);
Value       = zeros(nNetlist, 1);
node1       = cell(nNetlist, 1);
node2       = cell(nNetlist, 1);
for ii = 1:nNetlist
    cellNetlist{ii} = regexp(strNetlist{ii}, ' +', 'split');
    iType(ii)   = funSimType2iType(cellNetlist{ii}{2});
    cellName{ii}= cellNetlist{ii}{1};
    Value(ii)   = str2double(cellNetlist{ii}{5});
    node1{ii}   = cellNetlist{ii}{3};
    node2{ii}   = cellNetlist{ii}{4};
end
end