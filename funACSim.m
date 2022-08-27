%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
%--------------------------------------------------------------------------
function funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2)
freq = logspace(log10(f0), log10(f1), N);
Vo = zeros(1, N);
[a, bl]  = ismember('RL',cellName);
% [a, bs]   = ismember('RS',cellName);
[a2, b2] = ismember('iRL',MX);
RL  = Value(bl);
% RS  = Value(bs);
if RL == 0
    ib2 = b2;
else
    nRL = max(node1(bl), node2(bl));
    ib2 = nRL;
end
for ii=1:N
    f = freq(ii);
    s = 1i*2*pi*f;
    V = (MM.*s + MN)\MV';
    Vo(ii) = V(ib2);
end
dBVo = 20*log10(abs(Vo));
AgVo = angle(Vo)*180/pi;
uWVo = unwrap(AgVo, 179);
% uWVo = AgVo;
% toc;
% figure(2)
semilogx(axMag, freq, dBVo, '-r', 'LineWidth', 2);
hold(axMag, 'off');
grid(axMag, 'on');
xlabel(axMag, 'Freq/Hz');
if RL == 0
    ylabel(axMag, 'I_o Mag/dB');
else
    ylabel(axMag, 'V_o Mag/dB');
end
title(axMag, 'FrequencyResponse');
xlim(axMag, [min(freq),max(freq)]);
ylim(axMag, [-80,0]);

semilogx(axPhase, freq, uWVo, '-r', 'LineWidth', 2);
xlim(axPhase, [min(freq),max(freq)]);
ylim(axPhase, 'auto');
grid(axPhase, 'on');
xlabel(axPhase, 'Freq/Hz');
if RL == 0
    ylabel(axPhase, 'I_o Angle/deg');
else
    ylabel(axPhase, 'V_o Angle/deg');
end
title(axPhase, 'Phase VS. Freq');
end