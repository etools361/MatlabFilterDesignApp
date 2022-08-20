%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 找出主路径
%--------------------------------------------------------------------------
function [cellBranch2GND, iBranch2GND, Bn, iBn, e] = funFind2GNDPath(h, e, PathTemp, n)
% 标记主链路节点
NodeAll = zeros(1, n);
NodeAll(PathTemp) = 1;
% 找出对地枝节
cellBranch2GND = [];
iBranch2GND    = 0;
% 以ns为起点，以1节点（GND）为终点的枝节搜索路径,路径存储在cellBranch2GND中
ne = 1;
ins = 0;
HistPath   = [];
HistDevice = [];
% 深搜找出对地枝节
for ns = PathTemp
    ins = ins + 1;
    [NodeAll, cellBranch2GND, iBranch2GND, HistPath, HistDevice, e] = funDFSFindBranch(h, e, ns, ne, NodeAll, cellBranch2GND, iBranch2GND, HistPath, PathTemp(ins), HistDevice);
end
Bn = [];
for cP = 1:iBranch2GND
    Bn = [Bn, cellBranch2GND{cP}{3}];
end
Bn = [Bn, PathTemp];
Bn = unique(Bn);
iBn = zeros(1, n);
for ii=Bn
    iBn(ii) = 1;
end
end