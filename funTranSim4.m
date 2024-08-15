%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
% 求解方程： MN*X + MM*dX/dt = MV, dX/dt = (MV-MN*X)*MM^(-1)
%--------------------------------------------------------------------------
function [t, Vo]=funTranSim4(f1, Tmax, Nmax, MX, MM, MN, MV, retNode)
h=Tmax/Nmax;
t=linspace(0,Tmax,Nmax);   % 构造时间序列
E=1+square(2*pi*f1*t);% 构造脉冲激励信号
X = [];
X(:,1)=zeros(size(MX));
MNT  = inv(MM/h+MN)*h;
II   = eye(size(MM))*h;
MMTM = MNT*MM/h-II;
[a,  b]  = ismember('iV0',MX);
[a2, b2] = ismember(sprintf('V%dn', retNode),MX);
for n=2:Nmax
    MV(b)  = E(n-1);
%     X(:,n) = X(:,n-1) + 1/h*(MNT*MV' + MMTM*X(:,n-1));
    fK1 = (MNT*MV' + MMTM*X(:,n-1));
    ykp1 = X(:,n-1) + 1/h*fK1;
%     MV(b)  = E(n);
    fK2 = (MNT*MV' + MMTM*ykp1);
    X(:,n) = X(:,n-1) + 1/h*(fK1 + fK2)/2;
%     K1 = MNT*MV' + MMTM*X(:,n-1);
%     K2 = MNT*MV' + MMTM*(X(:,n-1)+1/h*K1/2);
%     K3 = MNT*MV' + MMTM*(X(:,n-1)+1/h*K2/2);
%     K4 = MNT*MV' + MMTM*(X(:,n-1)+1/h*K3);
%     X(:,n) = X(:,n-1) + 1/h/6*(K1+2*K2+2*K3+K4);
end
Vo = X(b2,:);
% plot(t, E, '-b', 'LineWidth', 2);
% plot(t, X(b2,:), '-r', 'LineWidth', 1);
% % hold('on');
% % hold('off');
% xlim([min(t),max(t)]);
% grid('on');
% xlabel('Time/s');
% title('V_o VS. t');
% legend({'Vi', 'Vo'}, 'location', 'northeast');

end

