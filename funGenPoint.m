%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成点
%--------------------------------------------------------------------------
function [x, y]=funGenPoint(axPlot, x, y, r)
hold(axPlot, 'on');
if r == 0
    plot(axPlot, x, y, '-ok', 'LineWidth', 1, 'MarkerFaceColor','k');
else
    plot(axPlot, y, x, '-ok', 'LineWidth', 1, 'MarkerFaceColor','k');
end
hold(axPlot, 'off');