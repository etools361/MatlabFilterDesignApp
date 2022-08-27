%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成点
%--------------------------------------------------------------------------
function [x, y]=funGenPoint(axPlot, x, y, r)
hold(axPlot, 'on');
if r == 0
    plot(axPlot, x, y, '-ok', 'LineWidth', 3);
    plot(axPlot, x, y, '-*k', 'LineWidth', 2);
else
    plot(axPlot, y, x, '-ok', 'LineWidth', 3);
    plot(axPlot, y, x, '-*k', 'LineWidth', 2);
end
hold(axPlot, 'off');