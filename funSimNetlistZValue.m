%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 重新转换节点标号
%--------------------------------------------------------------------------
function [ZValue] = funSimNetlistZValue(iType, Value, w)
    if iType == 0 % V
        ZValue = 0;
    elseif iType == 1 % I
        ZValue = inf;
    elseif iType == 2 % R
        ZValue = Value;
    elseif iType == 3 % L
        ZValue = w*Value;
    elseif iType == 4 % C
        ZValue = 1/(w*Value);
    end
end