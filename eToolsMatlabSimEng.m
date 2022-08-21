%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 通用无源网络仿真引擎
%--------------------------------------------------------------------------
% netlist, PI型，两端接载
tic;
fType = 'LPF';
n     = 11;
Rs    = 1;
Rl    = 1;
fp    = 1;
fs    = 2;
Ap    = -3;
As    = Ap-15;
TeeEn = 1;% TeeEn=0:PI, TeeEn=1:Tee
[strNetlist] = funSynthesisButterworthFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As);
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 2 3 0';
'C1 C 3 0 1';
'L2 L 3 4 2';
'C3 C 4 0 1';
'R2 R 4 5 0';
'RL R 5 0 1'
};
% T型，一端接载
strNetlist1 = {
'V0 V 1 GND 1';
'RS R 1 7 1';
'CS C 7 6 1';
'RS2 R 1 8 1';
'LS L 8 6 1';
'RP R 5 6 10';
'C3 C 3 x 1.333';
'C4 C 3 y 0.2';
'R5 R x y 1.1';
'R6 R y z 1.2';
'R7 R z GND 1.4';
'L2 L 6 3 1.5';
'L4 L 3 5 0.5';
'L5 L z 7 0.5';
'RL R 5 GND 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 2 3 0';
'C1 C 3 8 1';
'R3 R 8 9 1';
'L5 L 9 0 1';
'C4 C 8 0 0.2';
'L2 L 3 4 2';
'C3 C 4 0 1';
'R2 R 4 5 0.1';
'RL R 5 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 0 0.314631';
'C3 C 2 3 0.02996';
'L3 L 2 3 0.155487';
'C4 C 3 0 0.396551';
'C5 C 4 3 0.082683';
'L5 L 4 3 0.126308';
'C6 C 4 0 0.314631';
'RL R 4 0 1'
};
% Butterworth BPF
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 0 0.491582';
'L2 L 2 0 0.052049';
'C3 C 2 3 0.019881';
'L3 L 3 4 1.287';
'C4 C 4 0 1.591';
'L4 L 4 0 0.016084';
'C5 C 4 5 0.019881';
'L5 L 5 6 1.287';
'C6 C 6 0 0.491582';
'L6 L 6 0 0.052049';
'RL R 6 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 3 0.1';
'R2 R 2 3 0.2';
'C3 C 4 3 1';
'R3 R 4 3 0.2';
'C4 C 4 5 10';
'R4 R 4 5 0.2';
'C5 C 6 5 100';
'R5 R 6 5 0.2';
'RL R 6 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 3 2 1';
'R2 R 2 0 1';
'R3 R 1 3 10';
'RL R 3 0 1'
};
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


