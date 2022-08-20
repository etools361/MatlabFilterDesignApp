%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 找出主路径
%--------------------------------------------------------------------------
function [DistanceTemp, PathTemp, DeviceMain, DeviceEdgeTemp, e] = funFindMainPath(h, e, sb0, db0)
[DistanceTemp, PathTemp, DeviceMain, DeviceEdgeTemp, e] = funDijkstra(h, e, sb0, db0);% 找出主链路
end