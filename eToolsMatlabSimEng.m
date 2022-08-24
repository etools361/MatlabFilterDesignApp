%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 通用无源网络仿真引擎
%--------------------------------------------------------------------------
% netlist, PI型，两端接载
tic;
fType  = 'Butterworth';
fShape = 'HPF';
n     = 11;
Rs    = 1;
Rl    = 1;
fp    = 1;
fs    = 2;
Ap    = -3;
As    = Ap-15;
bw    = [];
TeeEn = 1;% TeeEn=0:PI, TeeEn=1:Tee
% 滤波器综合
[strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape);
% 滤波器仿真
% split netlist
[iType, Value, cellNode1, CellNode2, cellName] = funSimNetlist2Array(strNetlist);
% netlist standard
[node1, node2] = funSimNetlistRenode(cellNode1, CellNode2);
% netlist analysis
[maxNode, nL, nI, nV, nR0] = funSimNetlistAna(iType, Value, node1, node2);
% 标记给定器件的非GND节点用于获取结果位置
strDevice = 'RL'; % 显示结果的器件
[retNode] = funGetDeviceNode(cellName, node1, node2, strDevice);
% generate matrix, (MM*s+MN)*MX = MV
[MM, MN, MV, MX] = funSimNetlist2Matrix(iType, Value, node1, node2, maxNode, nL, nI, nV, nR0, cellName);

%% -----------------Gen Sch-------------------------
figure(1);
axPlot = gca;
f0     = 1e-3/2/pi;
w0     = 2*pi*f0;
[img] = funGenSchematic2(axPlot, iType, Value, cellNode1, CellNode2, cellName, w0);
% tic;
% ylimValue = ylim(axPlot);
% ylim(axPlot, [-ceil(-ylimValue(1)), ceil(ylimValue(2))]);
% set(axPlot,'ytick',[]);
% set(axPlot,'xtick',[]);
% set(axPlot,'Box','off');
set(gcf,'color',[1,1,1]);
% set(gca,'looseInset',[0 0 0 0]);
% axis off;
% set(gca,'LooseInset',get(gca,'TightInset'));
% set(gca,'PlotBoxAspectRatio', [8 2 1])
% figure(2);
%% -----------------AC-------------------------
f0 = 1e-2;
f1 = 1e1;
% f0 = 1e4;
% f1 = 1e6;
N     = 500;
h2    = figure(2);
axMag = axes(h2);
h3    = figure(3);
axPhase = axes(h3);
funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2);

%% -----------------Tran-------------------------
f1     = 2e-2;
% f1=1;
% Tmax=30/f1;
Tmax   = 2/f1;
Nmax   = 10000; % 
h4     = figure(4);
axTran = axes(h4);
oldXTick = findall(axTran,'Tag', 'pseudoXTickLabels'); 
delete(oldXTick)
funTranSim(axTran, f1, Tmax, Nmax, MX, MM, MN, MV, retNode, cellName, Value)
toc;


