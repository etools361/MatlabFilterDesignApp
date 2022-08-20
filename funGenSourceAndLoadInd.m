%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成源和负载索引号
%--------------------------------------------------------------------------
function [sb0, db0, eSource, eLoad] = funGenSourceAndLoadInd(e, t, cellName)
nt = cellfun(@(x)length(x)/2, t);
% 确定其起始位置和终止位置
% 0:V,1:I,2:R,3:L,4:C
if nt(1) ~= 0
    eSource = t{1}(1);
elseif nt(2) ~= 0
    eSource = t{2}(1);
else
    eSource = 1;
end
[a, b] = ismember('RL', cellName);
if a
    eLoad = b*2-1;
else
    eLoad   = t{3};% is R and R to GND
    for ii=1:nt(3)
        if eLoad(ii) == 0
            eLoad = eLoad(ii);
            break;
        end
    end
    if isempty(eLoad)
        eLoad = idx;
    end
end
sb0 = e(eSource).node1;
db0 = e(eLoad).node1;
end