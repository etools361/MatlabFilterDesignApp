%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 绘制节点和线段
%--------------------------------------------------------------------------
function [AxisOfSchNode, x, y] = funDrawNodeAndLine(axPlot, h, x, y, PointLength, AxisOfSchNode, r, m, ii, cNode)
nodeT = h{cNode};
AxisOfSchNode{cNode} = {1, [x, y]};% 标记
if ii == 0
    mPoint = length(nodeT) - 2;
else
    if ii == 1 || ii == m
        ExtMainPathLength(ii) = length(nodeT) - 1;
    else
        ExtMainPathLength(ii) = length(nodeT) - 2;
    end
    mPoint = ExtMainPathLength(ii);
end
for jj=1:mPoint
    if jj~=1
        if r==0
            x = x + PointLength;
            funGenLine(axPlot, [x-PointLength, x], [y, y], r);
        else
            y = y - PointLength;
            funGenLine(axPlot, [y-PointLength, y], [x, x], r);
        end
        AxisOfSchNode{cNode} = [AxisOfSchNode{cNode}, [x, y]];
    end
    if ((jj ~= 1) || (ii ~= 1)) && ((jj ~= mPoint) || (ii ~= m))
        funGenPoint(axPlot, x, y, 0);
    end
end
end