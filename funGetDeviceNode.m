%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 获取给定器件的节点标号，除掉GND节点
%--------------------------------------------------------------------------
function [retNode] = funGetDeviceNode(cellName, node1, node2, strDevice)
[a, b] = ismember(strDevice, cellName);
if a
    retNode1 = node1(b);
    retNode2 = node2(b);
    if retNode1 ~= 0
        retNode = retNode1;
    elseif retNode2 ~= 0
        retNode = retNode2;
    else
        retNode = 1;
        warning(sprintf('Device "%s" Error!', strDevice));
    end
else
    retNode = 1;
    warning(sprintf('Device "%s" Error!', strDevice));
end
end