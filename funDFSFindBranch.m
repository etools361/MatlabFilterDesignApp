%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-15(yyyy-mm-dd)
% 深度搜索枝节
%--------------------------------------------------------------------------
function [NodeAll, cellBranch2GND, iBranch2GND, HistPath, HistDevice, e] = funDFSFindBranch(h, e, ns, ne, NodeAll, cellBranch2GND, iBranch2GND, HistPath, ins, HistDevice)
cur_all_node = h{ns};
for cNode = cur_all_node
    next_node = e(cNode).node2;
    if NodeAll(next_node) == 0
        HistPath(ns)   = next_node;
        HistDevice{ns} = e(cNode);
        if ne == next_node && HistPath(ins)% 找到终点节点
            iBranch2GND = iBranch2GND + 1;
            cellBranch2GND{iBranch2GND}{1} = ins;
            PathBranch = ins;
            PathBranchEdge = [];
            for kk=2:length(HistPath)
                PathBranch(kk) = HistPath(PathBranch(kk-1));
                edgeTemp = HistDevice{PathBranch(kk-1)}.edge;
                e(edgeTemp).isMarked = 1;% 标记器件
                e(bitxor(edgeTemp-1, 1) + 1).isMarked = 1;% 标记器件
                PathBranchEdge{kk-1} = HistDevice{PathBranch(kk-1)};
                if PathBranch(kk) == ne
                    isFind = 1;
                    break;
                end
            end
            cellBranch2GND{iBranch2GND}{2} = HistPath;
            cellBranch2GND{iBranch2GND}{3} = PathBranch;
            cellBranch2GND{iBranch2GND}{4} = PathBranchEdge;
            HistPath = [];
        else
            NodeAll(next_node) = 1;
            ns0 = next_node;
            [NodeAll, cellBranch2GND, iBranch2GND, HistPath, HistDevice, e] = funDFSFindBranch(h, e, ns0, ne, NodeAll, cellBranch2GND, iBranch2GND, HistPath, ins, HistDevice);
        end
    end
end
end