%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-05-31(yyyy-mm-dd)
% 多项式代入计算
%--------------------------------------------------------------------------
% 将 p1代入到p2中
function result = substitute_poly(p2, p1)
% 初始化结果的系数为零
result = p2(end);
p1_power = 1; 

% 对 p2 的每一项进行处理
for i = length(p2)-1:-1:1
    % 计算 p1 的 i-1 次幂
    p1_power = conv(p1_power, p1);
    
    % 将 p1 的 i-1 次幂乘以 p2 的第 i 项，然后加到结果中
    result = polyadd(result, p2(i) * p1_power);
end
end
