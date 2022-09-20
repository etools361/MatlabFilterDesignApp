%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 由网表生成原理图
% 如何由网表生成原理图？
% 前提，需要给定起点(源)和终点(负载),其他电路均为待定电路网络，也即将整个电路理解为一个二端口网络。
%       0节点为接地点，不用考虑。
% 使用图论方法，建立从起点到终点的最长树，再在这棵树上寻每个枝节，可能存在主树上长出一个连通的小树。
% 数据结构：使用链式前向星无向图结构
%--------------------------------------------------------------------------
function [img] = funGenSchematic2(axPlot, iType, Value, cellNode1, CellNode2, cellName, w0)
% netlist standard
% 生成无向图
[h, t, e, idx, nType, n] = funCreateUndirectedGraphs(iType, Value, cellNode1, CellNode2, cellName);
% 生成源和负载
[sb0, db0, eSource, eLoad] = funGenSourceAndLoadInd(e, t, cellName);
% 添加w0处的权重
[e] = funAddWeight(e, iType, abs(Value), w0, nType);% 避免负值影响结果
% 找出主链路,主链路器件为DeviceMain, 节点为PathTemp
[DistanceTemp, PathTemp, DeviceMain, DeviceEdgeTemp, e] = funFindMainPath(h, e, sb0, db0);
% 找出对地路径
[cellBranch2GND, iBranch2GND, Bn, iBn, e]      = funFind2GNDPath(h, e, PathTemp, n);
% 找出带节点的桥
[iBn, cellLoop, cellDevice0, iLoop, Looping, e] = funFindBridgeWiNode(h, e, Bn, iBn);
% 剩下未被标记的即为无节点桥。无桥节点需要预占位Cross Point，不然构建的原理图不对称
[mB, cellWoBridge] = funFindBridgeWoNode(h, e, PathTemp, DeviceEdgeTemp, cellBranch2GND, cellLoop, cellDevice0);

% 策略2
% 寻找从指定两个节点的最长路径，使用图的深度搜索，研究了一下，貌似这个问题是NP-Hard问题！参考哈密顿路径。
% [histNode] = funGetLongestEdge(h, t, e, eSource, eLoad);

% 原理图绘制
% 计算主链路坐标
% 计算主链路节点扩展值
% 计算出每个节点的坐标，并记录，用于构建原理图
CrossInfo = zeros(1, n); % 记录节点类型, 1:主链路节点,2:对地枝节点,3:桥枝节点
CrossInfo(PathTemp) = 1;
m = length(PathTemp);
ExtMainPathLength = zeros(1, m);
for ii=1:m
    nodeT = h{PathTemp(ii)};
    if ii == 1 || ii == m
        ExtMainPathLength(ii) = length(nodeT)-1;
    else
        ExtMainPathLength(ii) = length(nodeT)-2;
    end
end
axis(axPlot, 'auto');
plot(axPlot, 0,0);
AxisOfSchNode = cell(1, n); % 标记所绘节点坐标
m = length(PathTemp);
x = 0;
y = 0;
PointLength  = 0.8;
DeviceLength = 1;
r = 0;
MarkedEndNode = 0;
% 绘制主链路
for ii=1:m
    cNode = PathTemp(ii);
    [AxisOfSchNode, x, y] = funDrawNodeAndLine(axPlot, h, x, y, PointLength, AxisOfSchNode, r, m, ii, cNode);
    % 占位
    if MarkedEndNode % 上一个器件的最后一个点是否占位，没有，则需要占第一个点
        AxisOfSchNode{cNode}{1} = AxisOfSchNode{cNode}{1} + 1;
        MarkedEndNode = 0;
    end
    % draw wo node bridge
    if cellWoBridge{cNode}{1}
        for jj=1:cellWoBridge{cNode}{1}
            eDevice = cellWoBridge{cNode}{jj+1};
            iTypeD = eDevice.iType;
            ValueD = eDevice.Value;
            funGenLine(axPlot, [y, y+jj], [x, x], ~r);
            funGenLine(axPlot, [y, y+jj], [x+1, x+1], ~r);
            funDrawDevice(axPlot, iTypeD, ValueD, x, y+jj, r);
            % 标记
            eTemp = eDevice.edge;
            e(eTemp).isMarked = 1;
            e(bitxor(eTemp-1,1)+1).isMarked = 1;
            % 需要在画的最后一个点上占位
            cTemp = AxisOfSchNode{cNode};
            Temp1 = cTemp{1};
            Temp2 = length(cTemp);
            if Temp2 > 2
                TempNew = [cTemp(1:Temp1),cTemp(Temp2), cTemp(Temp1+1:Temp2-1)];
                AxisOfSchNode{cNode} = TempNew;
            end
            AxisOfSchNode{cNode}{1} = AxisOfSchNode{cNode}{1} + 1;
            MarkedEndNode = eDevice.node2;% 记录最后一个点
        end
    end
    % draw main device
    if ii ~= m
        iTypeD = e(DeviceEdgeTemp(ii)).iType;
        ValueD = e(DeviceEdgeTemp(ii)).Value;
        funDrawDevice(axPlot, iTypeD, ValueD, x, y, r);
        x = x + DeviceLength;
    end
end
xLoad = x;
% 绘制枝节点，对地枝节,如遇到起始终止枝节则放置在两边
r = 1;
m = iBranch2GND;
for ii=1:m
    BTemp = cellBranch2GND{ii};
    cellDevice = BTemp{4};
    mDevice    = BTemp{3};
    mB = length(mDevice)-1;
    cellMark    = AxisOfSchNode{mDevice(1)};% 读取标记
    if cellDevice{1}.edge == eLoad && (length(cellMark)-cellMark{1})>1% 负载在最终位置，交换位置
        lastPoint = cellMark{end};
        cellMark{end}   = cellMark{cellMark{1}+1};
        cellMark{cellMark{1}+1} = lastPoint;
        AxisOfSchNode{mDevice(1)} = cellMark;
    end
    NodeAxis    = cellMark{cellMark{1}+1};
    x           = NodeAxis(1);
    y           = NodeAxis(2);
    for jj=1:mB % 将隐藏的节点坐标记录下来
        strctDevice = cellDevice{jj};
        iTypeD      = strctDevice.iType;
        ValueD      = strctDevice.Value;
        cNode       = mDevice(jj);
        if jj~=1
            y = y - DeviceLength;
            [AxisOfSchNode, x, y] = funDrawNodeAndLine(axPlot, h, x, y, PointLength, AxisOfSchNode, r, m, 0, cNode);
            CrossInfo(cNode) = 2;
        end
        funDrawDevice(axPlot, iTypeD, ValueD, y-DeviceLength, x, r);
    end
    funGenGND(axPlot, x, y-DeviceLength, r);
    AxisOfSchNode{mDevice(1)}{1} = AxisOfSchNode{mDevice(1)}{1} + 1;
end
% 绘制带节点桥
r = 0;
m = length(cellLoop);
for ii=1:m
    mLoop   = cellLoop{ii};
    mDevice = cellDevice0{ii};
    mD      = length(mDevice);
    % 获取坐标，取起始终止坐标，求取对应坐标位置
    AxisStart = AxisOfSchNode{mLoop(1)}{AxisOfSchNode{mLoop(1)}{1}+1};
    AxisStop  = AxisOfSchNode{mLoop(end)}{AxisOfSchNode{mLoop(end)}{1}+1};
    % 若不在同一水平线上，需要画线
    if AxisStart(2) ~= AxisStop(2)
        if abs(AxisStart(2)) > abs(AxisStop(2))
            y0 = AxisStop(2);
            y1 = AxisStart(2);
            x  = AxisStop(1);
        else
            y0 = AxisStart(2);
            y1 = AxisStop(2);
            x  = AxisStart(1);
        end
        funGenLine(axPlot, [y0, y1], [x, x], ~r); % 画垂直线
    else
        warning('error!');
        y1 = 0;
    end
    y = y1;
    dXMain = abs(AxisStart(1)-AxisStop(1));% 求主链路起始终止水平长度
    if dXMain>mD
        % 涉及到扩展桥长度
        dExpLength = (dXMain-mD)/mD;
    else
        dExpLength = 0;
        % 涉及到拉伸，后续考虑
        warning('需要拉伸!\n');
    end
    for jj=1:mD
        cDevice = mDevice{jj};
        iTypeD = cDevice.iType;
        ValueD = cDevice.Value;
        if jj ~= 1
            cNode = mLoop(jj);
            [AxisOfSchNode, x, y] = funDrawNodeAndLine(axPlot, h, x, y, PointLength, AxisOfSchNode, r, m, 0, cNode);
            CrossInfo(cNode) = 3;
        end
        if AxisStart(1) < AxisStop(1)
            x = x + DeviceLength;
        else
            x = x - DeviceLength;
        end
        funDrawDevice(axPlot, iTypeD, ValueD, x, y, r);
        if dExpLength
            if AxisStart(1) < AxisStop(1)
                x1 = x + dExpLength;
            else
                x1 = x - dExpLength;
            end
            funGenLine(axPlot, [x, x1], [y, y], r); % 画水平线
            x = x1;
        end
    end
    AxisOfSchNode{mLoop(1)}{1} = AxisOfSchNode{mLoop(1)}{1} + 1;
    AxisOfSchNode{mLoop(end)}{1} = AxisOfSchNode{mLoop(end)}{1} + 1;
end
% 绘制无节点桥
for ii=1:idx
    nDevice = e(ii);
    % 需要绘制的无节点桥
    if nDevice.isMarked == 0 && mod(ii, 2) == 1
        Pv = 0;
        node1 = nDevice.node1;
        node2 = nDevice.node2;
        if node2 == 1 % 处理一个节点为GND的特殊情况
            node1 = nDevice.node2;
            node2 = nDevice.node1;
        end
        N2 = AxisOfSchNode{node2};
        if node1 == 1
            N1 = N2;
            N1{2}(2) = N2{2}(2)-1;
        else
            N1 = AxisOfSchNode{node1};
        end
        P1 = N1{N1{1}+1};
        P2 = N2{N2{1}+1};
        if P1(2) ~= P2(2)
            if abs(P1(2)) < abs(P2(2))
                y0 = P1(2);
                y1 = P2(2);
                x1 = P1(1);
                x2 = P2(1);
            else
                y0 = P2(2);
                y1 = P1(2);
                x1 = P2(1);
                x2 = P1(1);
            end
            if x1 == x2 % 属于垂直桥
                Pv = 1;
                if CrossInfo(node2) == 2 % 对地枝节上的才需要横向扩展
                    x = x1+0.5;
                    funGenLine(axPlot, [x1, x], [y0, y0], r); % 画水平线
                else
                    x = x1;
                end
                if node1 ~= 1 % 不是到GND的，需要画横线
                    funGenLine(axPlot, [x1, x], [y1, y1], r); % 画水平线
                end
            else % 属于水平桥
                if CrossInfo(node1) == 3 || CrossInfo(node2) == 3
                    funGenLine(axPlot, [y1, y1-1], [x2, x2], ~r); % 画垂直线
                    y1 = y1 - 1;
                end
                funGenLine(axPlot, [y0, y1], [x1, x1], ~r); % 画垂直线
            end
        else
            % 往上或下一个单位，具体需要由节点性质决定，若CrossInfo=1，需要往上，CrossInfo=3需要往下
            y0 = P1(2);
            if CrossInfo(node1) == 1
                y1 = y0 + 1;
            else
                y1 = y0 - 1;
            end
            x0 = P1(1);
            x1 = P2(1);
            funGenLine(axPlot, [y0, y1], [x0, x0], ~r); % 画垂直线
            funGenLine(axPlot, [y0, y1], [x1, x1], ~r); % 画垂直线
        end
        if Pv
            cDevice = nDevice;
            iTypeD = cDevice.iType;
            ValueD = cDevice.Value;
            funDrawDevice(axPlot, iTypeD, ValueD, y1, x, ~r);
            if node1 == 1 % 到GND的
                funGenGND(axPlot, x, y1, ~r);
            else
            end
        else
            dXMain = abs(P1(1)-P2(1));% 求主链路起始终止水平长度
            if dXMain>1
                % 涉及到扩展桥长度
                dExpLength = (dXMain-1)/1;
            else
                dExpLength = 0;
                % 涉及到拉伸，后续考虑
            end
            cDevice = nDevice;
            iTypeD = cDevice.iType;
            ValueD = cDevice.Value;
    %         funGenLine([x, x], [y0, y1], r); % 画垂直线
            y  = y1;
            if P1(1) < P2(1)
                x1 = P1(1);
                x2 = P2(1);
            else
                x1 = P2(1);
                x2 = P1(1);
            end
            funGenLine(axPlot, [y, y], [x1, x1+dExpLength/2], ~r); % 画水平直线
            funGenLine(axPlot, [y, y], [x2-dExpLength/2, x2], ~r); % 画水平直线
            funDrawDevice(axPlot, iTypeD, ValueD, x1+dExpLength/2, y, r);
%             fprintf('DeviceId = %d\n', ii);
            if node1 ~= 1
                AxisOfSchNode{node1}{1} = AxisOfSchNode{node1}{1} + 1;
            end
            AxisOfSchNode{node2}{1} = AxisOfSchNode{node2}{1} + 1;
        end
    end
end
axis(axPlot, 'equal');
hold(axPlot, 'off');

x0 = xLoad;
y0 = 0;
[a, b]   = ismember('RL',cellName);
RL  = Value(b);
if RL == 0
    text(axPlot,x0+0.0, y0+0.15, 'I_o', 'FontSize',12, 'FontWeight', 'bold');
else
    text(axPlot,x0+0.0, y0+0.15, 'V_o', 'FontSize',12, 'FontWeight', 'bold');
end
img = [];
set(axPlot,'ytick',[]);
set(axPlot,'xtick',[]);
set(axPlot,'Box','off');
axis(axPlot, 'off');
ylimValue = ylim(axPlot);
ylim(axPlot, [-ceil(-ylimValue(1)*5)/5, ceil(ylimValue(2)*5)/5]);


function [x1, y1] = funDrawDevice(axPlot, iType, Value, x2, y2, r)
% 0:V,1:I,2:R,3:L,4:C
switch iType
    case 0 % V
        [x1, y1] = funGenV(axPlot, x2, y2, r);
        funGenText(axPlot, x2+1.1, y2-0, r, Value, 'V');
        text(axPlot, x2+0.8, y2+0.18, 'V_i', 'FontSize',12, 'FontWeight', 'bold');
        y0 = y1;
    case 1 % I
        [x1, y1] = funGenI(axPlot, x2, y2, r);
        funGenText(axPlot, x2+0.8, y2-0.55, r, Value, 'A');
        text(axPlot, x2+1, y2+0.18, 'I_i', 'FontSize',12, 'FontWeight', 'bold');
        y0 = y1;
    case 2 % R
        if Value == 0
            [x1, y1] = funGenLine(axPlot, [x2,x2+1], [y2,y2], r);
        elseif Value > 1e24
            [x1, y1] = funGenROpen(axPlot, x2, y2, r);
        else
            [x1, y1] = funGenR(axPlot, x2, y2, r);
            funGenText(axPlot, x1, y1, r, Value, '\Omega');
        end
    case 3 % L
        [x1, y1] = funGenL(axPlot, x2, y2, r);
        funGenText(axPlot, x1, y1, r, Value, 'H');
    case 4 % C
        [x1, y1] = funGenC(axPlot, x2, y2, r);
        funGenText(axPlot, x1, y1, r, Value, 'F');
end


function [x, y]=funGenText(axPlot, x, y, r, Value, strUnit)
hold(axPlot, 'on');
if r == 0
%     text(x-0.1, y+0.15, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',12, 'FontWeight', 'bold');
    t = text(axPlot, x-0.7, y+0.35, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',9, 'FontWeight', 'bold');
else
%     text(x+0.1, y-0.9, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',12, 'FontWeight', 'bold');
    t = text(axPlot, x+0.1, y-0.95, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',9, 'FontWeight', 'bold');
end
% t.Interactions = editInteraction;
hold(axPlot, 'off');


function [x0, y0]=funGenR(axPlot, x, y, r)
rb = 0.2;
ll = 1.0;
rh = 0.25;
rv = ll-2*rb;
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], [y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], '-k', 'LineWidth', 2);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
    plot(axPlot, [y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], [x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], '-k', 'LineWidth', 2);
    x0 = y;
    y0 = x+ll;
end
hold(axPlot, 'off');

function [x0, y0]=funGenROpen(axPlot, x, y, r)
rb = 0.2;
ll = 1.0;
rh = 0.25;
rv = ll-2*rb;
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
%     plot([x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], [y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
%     plot([y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], [x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x+ll;
end
hold(axPlot, 'off');

function [x0, y0]=funGenC(axPlot, x, y, r)
rb = 0.45;
ll = 1.0;
rh = 0.5;
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+rb, x+rb], [y+rh/2, y-rh/2], '-k', 'LineWidth', 2);
    plot(axPlot, [x+ll-rb, x+ll-rb], [y+rh/2, y-rh/2], '-k', 'LineWidth', 2);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
    plot(axPlot, [y+rh/2, y-rh/2], [x+rb, x+rb], '-k', 'LineWidth', 2);
    plot(axPlot, [y+rh/2, y-rh/2], [x+ll-rb, x+ll-rb], '-k', 'LineWidth', 2);
    x0 = y;
    y0 = x+ll;
end
hold(axPlot, 'off');

function [x0, y0]=funGenL(axPlot, x, y, r)
rb = 0.18;
ll = 1.0;
r0 = ll-rb*2;
lx = linspace(x+rb, x+ll-rb, 49);
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, lx, y+0.2.*abs(sin((lx-rb-x)./r0.*4.*pi)).^0.5-0.005, '-k', 'LineWidth', 2);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
    plot(axPlot, y+0.2.*abs(sin((lx-rb-x)./r0.*4.*pi)).^0.5-0.005, lx, '-k', 'LineWidth', 2);
    x0 = y;
    y0 = x+ll;
end
hold(axPlot, 'off');

function [x0, y0]=funGenGND(axPlot, x, y, r)
rb = 0.18;
dy = 0.1;
dx0 = 0.4;
dx1 = 0.25;
dx2 = 0.1;
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [y, y-rb], [x, x], '-k', 'LineWidth', 1);
    plot(axPlot, [y-rb, y-rb], [x-dx0/2, x+dx0/2], '-k', 'LineWidth', 2);
    plot(axPlot, [y-rb-dy, y-rb-dy], [x-dx1/2, x+dx1/2], '-k', 'LineWidth', 2);
    plot(axPlot, [y-rb-dy*2, y-rb-dy*2], [x-dx2/2, x+dx2/2], '-k', 'LineWidth', 2);
    x0 = y;
    y0 = x;
else
    plot(axPlot, [x, x], [y, y-rb], '-k', 'LineWidth', 1);
    plot(axPlot, [x-dx0/2, x+dx0/2], [y-rb, y-rb], '-k', 'LineWidth', 2);
    plot(axPlot, [x-dx1/2, x+dx1/2], [y-rb-dy, y-rb-dy], '-k', 'LineWidth', 2);
    plot(axPlot, [x-dx2/2, x+dx2/2], [y-rb-dy*2, y-rb-dy*2], '-k', 'LineWidth', 2);
    x0 = x;
    y0 = y;
end
hold(axPlot, 'off');

function [x0, y0]=funGenV(axPlot, x, y, r)
rb = 0.18;
ll = 1.0;
d0 = ll-rb*2;
d1 = 0.2;
d2 = 0.15;
lx = linspace(0, 2*pi, 41);
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, d0/2.*cos(lx)+x+ll/2, d0/2.*sin(lx)+y, '-k', 'LineWidth', 2);
    plot(axPlot, [x+ll-rb-d1+d2, x+ll-rb-d1], [y, y], '-k', 'LineWidth', 2);
    plot(axPlot, [x+ll-rb-d1+d2/2, x+ll-rb-d1+d2/2], [y-d2/2, y+d2/2], '-k', 'LineWidth', 2);
    plot(axPlot, [x+rb+d1-d2/2, x+rb+d1-d2/2], [y-d2/2, y+d2/2], '-k', 'LineWidth', 2);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
    plot(axPlot, d0/2.*sin(lx)+y, d0/2.*cos(lx)+x+ll/2, '-k', 'LineWidth', 2);
    plot(axPlot, [y, y], [x+ll-rb-d1+d2, x+ll-rb-d1], '-k', 'LineWidth', 2);
    plot(axPlot, [y-d2/2, y+d2/2], [x+ll-rb-d1+d2/2, x+ll-rb-d1+d2/2], '-k', 'LineWidth', 2);
    plot(axPlot, [y-d2/2, y+d2/2], [x+rb+d1-d2/2, x+rb+d1-d2/2], '-k', 'LineWidth', 2);
    x0 = x;
    y0 = y+ll;
end
hold(axPlot, 'off');

function [x0, y0]=funGenI(axPlot, x, y, r)
rb = 0.18;
ll = 1.0;
d0 = ll-rb*2;
d1 = 0.2;
d2 = 0.15;
lx = linspace(0, 2*pi, 41);
hold(axPlot, 'on');
if r == 0
    plot(axPlot, [x, x+rb], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, [x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 1);
    plot(axPlot, d0/2.*cos(lx)+x+ll/2, d0/2.*sin(lx)+y, '-k', 'LineWidth', 2);
    plot(axPlot, [x+ll-rb-d1+d2, x+rb+d1-d2], [y, y], '-k', 'LineWidth', 2);
    plot(axPlot, [x+ll-rb-d1, x+ll-rb-d1+d2, x+ll-rb-d1], [y-d2/2, y, y+d2/2], '-k', 'LineWidth', 2);
    x0 = x+ll;
    y0 = y;
else
    plot(axPlot, [y, y], [x, x+rb], '-k', 'LineWidth', 1);
    plot(axPlot, [y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 1);
    plot(axPlot, d0/2.*sin(lx)+y, d0/2.*cos(lx)+x+ll/2, '-k', 'LineWidth', 2);
    plot(axPlot, [y, y], [x+ll-rb-d1+d2, x+rb+d1-d2], '-k', 'LineWidth', 2);
    plot(axPlot, [y-d2/2, y, y+d2/2], [x+ll-rb-d1, x+ll-rb-d1+d2, x+ll-rb-d1], '-k', 'LineWidth', 2);
    x0 = x;
    y0 = y+ll;
end
hold(axPlot, 'off');

function [nodeLen, nodeNew, nodeChT, index, histNode] = funFindNodeNext(nodeChT, node1, node2, nodeLen, index, histNode, nodeS)
    itoNode = [find(nodeChT == node1),find(nodeChT == node2)];
    index = index + 1;
    nodeLen(itoNode) = index;
    nodeAll = [node1(itoNode),node2(itoNode)];
    histNode(index) = nodeChT;
    nodeNew = nodeAll(~ismember(nodeAll, [histNode,0]));

    nN = length(nodeNew);
    for jj=1:nN
        nodeChT = nodeNew(jj);
        if nodeChT ~= nodeS
            [nodeLen, nodeNew, nodeChT, index, histNode] = funFindNodeNext(nodeChT, node1, node2, nodeLen, index, histNode, nodeS);
        else
            histNode(index+1) = nodeChT;
            break;
        end
    end

function    [iNode, mNodeBefore, mNodeBack] = funFindNode(mT1, mT2, node1, node2)
findN1 = find(node1 == mT1);
findN2 = find(node2 == mT2);
if ~isempty(findN2)
    iNode1 = intersect(findN1, findN2);
else
    iNode1 = [];
end
findN1 = find(node1 == mT2);
findN2 = find(node2 == mT1);
if ~isempty(findN2)
    iNode2 = findN1(findN1 == findN2);
else
    iNode2 = [];
end
iNode = [iNode1, iNode2];
% 前面有多少个节点
NodeAll = [node1, node2];
ix = find(mT1 == NodeAll);
mNodeBefore = length(ix);
% 后面有多少个节点
ix = find(mT2 == NodeAll);
mNodeBack = length(ix);


