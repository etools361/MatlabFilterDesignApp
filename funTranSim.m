%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
% 求解方程： MN*X + MM*dX/dt = MV, dX/dt = (MV-MN*X)*MM^(-1)
%--------------------------------------------------------------------------
function funTranSim(axTran, f1, Tmax, Nmax, MX, MM, MN, MV, retNode, cellName, Value)
h = Tmax/Nmax;
t = linspace(0,Tmax,Nmax);   % 构造时间序列
X = [];
X(:,1)=zeros(size(MX));
MNT  = inv(MM/h+MN);
% II   = eye(size(MM))*h;
% MMTM = MNT*MM/h-II;
[a3, b3]  = ismember('RL', cellName);
[a4, b4]  = ismember('RS', cellName);
[a5, b5]  = ismember('V0', cellName);
I0 = 0;
if a5 == 0
    [a5, b5]  = ismember('I0', cellName);
    I0 = 1;
end
RL  = Value(b3);
RS  = Value(b4);
V0  = Value(b5);
if RL == 0 % current
    [a2, b2] = ismember('iRL',MX);
    E = RS*(1+square(2*pi*f1*t));% 构造脉冲激励信号
else % voltage
    [a2, b2] = ismember(sprintf('V%dn', retNode),MX);
    E = V0*(1+square(2*pi*f1*t));% 构造脉冲激励信号
end
[a,  b]  = ismember('iV0',MX);
if a==0
    [a,  b]  = ismember('iI0',MX);
end
for n=2:Nmax
    MV(b)  = E(n-1);
    X(:,n) = MNT*(MV' + 1/h*MM*X(:,n-1));
end
if RL == 0 % current
    plot(axTran, t, X(b2,:), '-r', 'LineWidth', 2);
    hold(axTran, 'off');
else
    if ~I0
        plot(axTran, t, E, '-b', 'LineWidth', 2);
        hold(axTran, 'on');
    end
    plot(axTran, t, X(b2,:), '-r', 'LineWidth', 2);
    hold(axTran, 'off');
end

xlim(axTran, [min(t),max(t)]);
ylim(axTran, 'auto');
grid(axTran, 'on');
xlabel(axTran, 'Time/s');
if RL == 0
    ylabel(axTran, 'I_o/A');
    title(axTran, 'I_o VS. t');
    legend(axTran, {'Io'}, 'location', 'northeast');
else
    ylabel(axTran, 'V_o/V');
    title(axTran, 'V_o VS. t');
    if ~I0
        legend(axTran, {'Vi', 'Vo'}, 'location', 'northeast');
    else
        legend(axTran, {'Vo'}, 'location', 'northeast');
    end
end

end

