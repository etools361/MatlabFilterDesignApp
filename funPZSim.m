%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
%--------------------------------------------------------------------------
function funPZSim(axPZ, axPZ3D, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, A0, Ap)
freq = linspace(-f1, f1, N);
[FX, FY] = meshgrid(freq, freq);
F        = FX+1i*FY;
Vo       = zeros(N, N);
[a,   b] = ismember('RL',cellName);
[a2, b2] = ismember('iRL',MX);
RL       = Value(b);
nRL      = max(node1(b), node2(b));
if RL == 0
    ib = b2;
else
    ib = nRL;
end
for jj=1:N
    for ii=1:N
        f = F(ii, jj);
        s = 2*pi*f;
        V = (MM*s + MN)\MV';
        Vo(ii, jj) = V(ib);
    end
end
Vo2 = Vo;
for jj=1:N
    for ii=1:N
        Vo2(ii, jj) = ((Vo(ii, jj)*Vo(ii, N+1-jj)));
    end
end
dBVo = 10*log10(abs(Vo2));
% dBVo(dBVo<-80) = NaN;
% fH = @(s) 1./(-s.^(2*n)+1);% 巴特沃斯滤波器
% H  = fH(F).*A0;
% Hi = 10*log10(abs(H));
% if Disp3DEn
axPZ1 = axPZ3D;
    surfc(axPZ1, FX, FY, dBVo, 'LineWidth', 0.1);
%     shading(axPZ1, 'interp');
%     lighting(axPZ1, 'flat');
%     light(axPZ1, 'position',[-1,0.5,1]*f1,'style','local','color','w')
    hold(axPZ1, 'on');
    % surf(axPZ, FX, FY, Hi);%ssss
    b  = freq;
    ibb = b>0;
    bb = b(ibb);
    plot3(axPZ1, bb.*0, bb, dBVo(ibb, fix(N/2)), '-*r', 'LineWidth', 2);
    hold(axPZ1, 'off');

    axis(axPZ1, 'tight');
    view(axPZ1, [1 1 1]);%侧视
    xlabel(axPZ1, 'real');
    ylabel(axPZ1, 'image');
    minZ = min(min(dBVo));
    if minZ>-80
        zlim(axPZ1, [minZ, A0+3]);
    else
        zlim(axPZ1, [-80, A0+3]);
    end
    xlim(axPZ1, [freq(1), freq(end)]);
    ylim(axPZ1, [freq(1), freq(end)]);
    grid(axPZ1, 'on');
    title(axPZ1, '3D PZ Plot');
% else
    contour(axPZ, FX, FY, dBVo', [minZ:(minZ-Ap),(A0+Ap):A0], 'LineWidth', 2);
%     hold on;
%     plot(axPZ, [1], [1], '*');
%     hold off;
    grid(axPZ,   'on');
    axis(axPZ,   'square');
    ylabel(axPZ, 'real');
    xlabel(axPZ, 'image');
    xlim(axPZ, [freq(1), freq(end)]);
    ylim(axPZ, [freq(1), freq(end)]);
    title(axPZ, 'PZ Plot');
% end

end