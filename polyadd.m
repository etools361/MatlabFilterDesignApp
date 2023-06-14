%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-05-31(yyyy-mm-dd)
% 多项式相加
%--------------------------------------------------------------------------
% 定义一个函数，用于计算两个多项式的和
function result = polyadd(p1, p2)
    % 如果 p1 的长度小于 p2，那么在 p1 前面补零
    if length(p1) < length(p2)
        p1 = [zeros(1, length(p2) - length(p1)) p1];
    % 否则，如果 p2 的长度小于 p1，那么在 p2 前面补零
    elseif length(p2) < length(p1)
        p2 = [zeros(1, length(p1) - length(p2)) p2];
    end
    
    % 计算 p1 和 p2 的和
    result = p1 + p2;
end