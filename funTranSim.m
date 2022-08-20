%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
% 求解方程： MN*X + MM*dX/dt = MV, dX/dt = (MV-MN*X)*MM^(-1)
%--------------------------------------------------------------------------
function funTranSim(axTran, f1, Tmax, Nmax, MX, MM, MN, MV, retNode, cellName, Value)
h=Tmax/Nmax;
t=linspace(0,Tmax,Nmax);   % 构造时间序列
E=1+square(2*pi*f1*t);% 构造脉冲激励信号
X = [];
X(:,1)=zeros(size(MX));
MNT  = inv(MM/h+MN);
% II   = eye(size(MM))*h;
% MMTM = MNT*MM/h-II;
[a,  b]  = ismember('iV0',MX);
[a2, b2] = ismember(sprintf('V%dn', retNode),MX);
for n=2:Nmax
    MV(b)  = E(n-1);
    X(:,n) = MNT*(MV' + 1/h*MM*X(:,n-1));
end
plot(axTran, t, E, '-b', 'LineWidth', 2);
hold(axTran, 'on');
plot(axTran, t, X(b2,:), '-r', 'LineWidth', 2);
hold(axTran, 'off');
xlim(axTran, [min(t),max(t)]);
ylim(axTran, 'auto');
grid(axTran, 'on');
xlabel(axTran, 'Time/s');
[a, b]   = ismember('RL', cellName);
RL  = Value(b);
if RL == 0
    ylabel(axTran, 'I_o/A');
else
    ylabel(axTran, 'V_o/V');
end
title(axTran, 'V_o VS. t');
legend(axTran, {'Vi', 'Vo'}, 'location', 'northeast');

end

