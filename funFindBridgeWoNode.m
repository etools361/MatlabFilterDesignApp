%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 找出不带节点的桥，并标记主链路上
% cellWoBridge：标记主链路上的桥(不带节点)个数
%--------------------------------------------------------------------------
function [mB, cellWoBridge] = funFindBridgeWoNode(h, e, PathTemp, DeviceEdgeTemp, cellBranch2GND, cellLoop, cellDevice0)
    n = length(e);
    WoBridgeMark = zeros(1, n);
    m = length(DeviceEdgeTemp);
    for ii=1:m
        mEdge = DeviceEdgeTemp(ii);
        WoBridgeMark(mEdge)                = 1;
        WoBridgeMark(bitxor(mEdge-1, 1)+1) = 1;
    end
    m = length(cellBranch2GND);
    for ii=1:m
        Branch2GND = cellBranch2GND{ii};
        n = length(Branch2GND{4});
        for jj=1:n
            eB = Branch2GND{4}{jj};
            mEdge = eB.edge;
            WoBridgeMark(mEdge)                = 1;
            WoBridgeMark(bitxor(mEdge-1, 1)+1) = 1;
        end
    end
    iBB = ~WoBridgeMark;
    cellWoBridge0 = e(iBB);
    nn = length(cellWoBridge0);
    cellWoBridge = cell(1, max(PathTemp));
    m = length(PathTemp);
    for ii=1:max(PathTemp)
        cellWoBridge{ii}{1} = 0;
    end
    mB = 0;
    for ii=1:m-1
        n1 = PathTemp(ii);
        n2 = PathTemp(ii+1);
        for jj=1:nn
            eTemp = cellWoBridge0(jj);
            if n1==eTemp.node1 && n2==eTemp.node2
                mB = mB + 1;
                cellWoBridge{n1}{1} = cellWoBridge{n1}{1} + 1;
                cellWoBridge{n1}{cellWoBridge{n1}{1}+1} = eTemp;
            end
        end
    end
end