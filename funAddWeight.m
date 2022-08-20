%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成源和负载索引号
%--------------------------------------------------------------------------
function [e] = funAddWeight(e, iType, Value, w0, nType)
for ii = 1:nType
    e(ii*2-1).ZValue = funSimNetlistZValue(iType(ii), Value(ii), w0);
    e(ii*2).ZValue   = funSimNetlistZValue(iType(ii), Value(ii), w0);
end
end