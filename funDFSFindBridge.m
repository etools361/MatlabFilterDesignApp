%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-15(yyyy-mm-dd)
% 深度搜索桥
% nn:节点判重
%--------------------------------------------------------------------------
function [iBn, cellLoop, cellDevice, iLoop, Looping, e] = funDFSFindBridge(h, e, ist, iBn, cellLoop, cellDevice, iLoop, Looping)
    next_node = h{ist};
    m_nn = length(next_node);
    i_m_nn = 0;
    for current_node = next_node
        nextN = e(current_node).node2;
        if iBn(nextN) == 0
            iBn(nextN) = 1;
            if Looping == 0
                Looping = 1;
                iLoop = iLoop + 1;
                cellLoop{iLoop} = ist;
                cellDevice{iLoop} = {};
            end
            % 标记
            e(current_node).isMarked = 1;
            e(bitxor(current_node-1, 1)+1).isMarked = 1;
            cellLoop{iLoop} = [cellLoop{iLoop}, nextN];
            cellDevice{iLoop} = [cellDevice{iLoop}, e(current_node)];
            % 迭代
            [iBn, cellLoop, cellDevice, iLoop, Looping, e] = funDFSFindBridge(h, e, nextN, iBn, cellLoop, cellDevice, iLoop, Looping);
        else
            i_m_nn = i_m_nn + 1;
        end
    end
    % 这个节点为死节点,那么构成一个环，即得到所求桥
    if i_m_nn == m_nn
        % 闭合节点
        if Looping == 1
            Looping = 0;
            nodeT = zeros(1,length(h));
            for ii=cellLoop{iLoop}
                nodeT(ii) = 1;
            end
            for ii=next_node
                if nodeT(e(ii).node2) == 0
                    e(ii).isMarked = 1;
                    e(bitxor(ii-1, 1)+1).isMarked = 1;
                    cellLoop{iLoop}   = [cellLoop{iLoop}, e(ii).node2];
                    cellDevice{iLoop} = [cellDevice{iLoop}, e(ii)];
                    break;
                end
            end
        end
    end
end