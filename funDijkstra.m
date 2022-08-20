%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-15(yyyy-mm-dd)
% Dijkstra算法
%--------------------------------------------------------------------------
function [DistanceTemp, PathTemp, DeviceMain, DeviceEdgeTemp, e] = funDijkstra(h, e, sb0, db0)
n             = length(h);
visited(1:0)  = 0;
distance(1:n) = inf;
distance(sb0) = 0;
visited(1)    = 1;% 排除GND节点
visited(sb0)  = 1;
u             = sb0;
parent(1:0)   = 0;
DeviceMain    = [];
DeviceEdge    = [];
for ii = 1:n-1
    id = h{u};
    for v = id
        if v~=1 && (e(v).ZValue+distance(u) < distance(e(v).node2))
            distance(e(v).node2) = e(v).ZValue+distance(u);
            parent(e(v).node2)   = u;
            DeviceMain{e(v).node2} = e(v);
            DeviceEdge(e(v).node2) = v;
        end
    end
    temp             = distance;
    temp(visited==1) = inf;
    [t,u]            = min(temp);
    visited(u)       = 1;
end
PathTemp   = [];
DeviceTemp = [];
DeviceEdgeTemp = [];
if parent(db0) ~= 0
    t        = db0;
    PathTemp = [db0];
    while t~=sb0
        P              = parent(t);
        PathTemp       = [P PathTemp];
        DeviceTemp     = [DeviceMain{t},DeviceTemp];
        e(DeviceEdge(t)).isMarked = 1;
        e(bitxor(DeviceEdge(t)-1, 1)+1).isMarked = 1;
        DeviceEdgeTemp = [DeviceEdge(t),DeviceEdgeTemp];
        t=P;
    end
    DeviceMain = DeviceTemp;
end
DistanceTemp = distance(db0);
end