%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成线段
%--------------------------------------------------------------------------
function [x0,y0]=funGenLine(axPlot, x, y, r)
hold(axPlot, 'on');
if r == 0
    plot(axPlot, x, y, '-k', 'LineWidth', 1);
    x0 = x(end);
    y0 = y(end);
else
    plot(axPlot, y, x, '-k', 'LineWidth', 1);
    x0 = y(end);
    y0 = x(end);
end
hold(axPlot, 'off');