%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 通用无源网络仿真引擎
%--------------------------------------------------------------------------
% netlist, PI型，两端接载
tic;
fType  = 'Gaussian';%'Bessel';%'Butterworth';
fShape = 'LPF';
n     = 3;
Rs    = 1;
Rl    = 1;
fp    = 1/2/pi;
fs    = 2;
Ap    = -3.01;
As    = -Ap-15;
bw    = [];
TeeEn = 1;% TeeEn=0:PI, TeeEn=1:Tee
SimSW = 3;
% 滤波器综合
[strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape);
% 滤波器仿真
% split netlist
[iType, Value, CellNode1, CellNode2, cellName] = funSimNetlist2Array(strNetlist);
% 实际器件替换
% [iType, Value] = funSynthesisFilterActual(iType, Value, 'E24', 'E24');
[strNetlist] = funSimArray2Netlist(iType, Value, CellNode1, CellNode2, cellName);
% strNetlist = {
% 'V0 V 9 0 1.000000e+00';
% 'RS R 1 9 1.000000e+00';
% 'C1 C 1 0 2.26';
% 'L2 L 1 2 1.51';
% 'C3 C 2 0 1.51';
% 'L4 L 2 3 2.26';
% 'RL R 3 0 1.000000e+00';
% };
% netlist standard
[node1, node2] = funSimNetlistRenode(CellNode1, CellNode2);
% netlist analysis
[maxNode, nL, nI, nV, nR0] = funSimNetlistAna(iType, Value, node1, node2);
% 标记给定器件的非GND节点用于获取结果位置
strDevice = 'RL'; % 显示结果的器件
[retNode] = funGetDeviceNode(cellName, node1, node2, strDevice);
% generate matrix, (MM*s+MN)*MX = MV
[MM, MN, MV, MX] = funSimNetlist2Matrix(iType, Value, node1, node2, maxNode, nL, nI, nV, nR0, cellName);

%% -----------------Gen Sch-------------------------
if bitand(SimSW, 1)
    figure(1);
    axPlot = gca;
    f0     = 1e-3/2/pi;
    w0     = 2*pi*f0;
    [img] = funGenSchematic2(axPlot, iType, Value, CellNode1, CellNode2, cellName, w0);
    set(gcf,'color',[1,1,1]);
end
%% -----------------AC-------------------------
if bitand(SimSW, 2)
    f0 = 1e-2;
    f1 = 1e0;
    % f0 = 1e4;
    % f1 = 1e6;
    N     = 500;
    h2    = figure(2);
    axMag = axes(h2);
    h3    = figure(3);
    axPhase = axes(h3);
    [IdealFreq, IdealMag, IdealPhase, P, Z] = funSimFilterIdeal(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape, f0, f1, N);
    funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase);
end
%% -----------------Tran-------------------------
if bitand(SimSW, 4)
    f1     = 2e-2;
    
    % f1=1;
    % Tmax=30/f1;
    Tmax   = 2/f1;
    Nmax   = 10000; % 
    h4     = figure(4);
    axTran = axes(h4);
    oldXTick = findall(axTran,'Tag', 'pseudoXTickLabels'); 
    delete(oldXTick)
    funTranSim(axTran, f1, Tmax, Nmax, MX, MM, MN, MV, retNode, cellName, Value);
end
%% -----------------PZ-------------------------
if bitand(SimSW, 8)
    f1       = 2;
    f0       = 1e-2;
    h5       = figure(5);
    axPZ     = axes(h5);
    h6       = figure(6);
    axPZ3D   = axes(h6);
    N        = 16*n;
    Disp3DEn = 1;
    A0       = 10*log10(1/4*sqrt(Rs/Rl));
    funPZSim(axPZ, axPZ3D, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, A0, Ap, P, Z);
    toc;
end

