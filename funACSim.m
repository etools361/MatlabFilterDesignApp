%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-18(yyyy-mm-dd)
% AC仿真
%--------------------------------------------------------------------------
function funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase, logscaleEn)
if logscaleEn
    freq = logspace(log10(f0), log10(f1), N);
else
    freq = linspace(f0, f1, N);
end
Vo = zeros(1, N);
Ii = zeros(1, N);
[a, bl]  = ismember('RL',cellName);
[a0, b0]  = ismember('RS',cellName);
% [a, bs]   = ismember('RS',cellName);
[a2, b2] = ismember('iRL',MX);
[a3, b3] = ismember('iV0',MX);
if b3 == 0
    [a3, b3] = ismember('iI0',MX);
end
RL  = Value(bl);
RS  = Value(b0);
ib3 = b3;
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
    Ii(ii) = V(ib3);
end
S11r = real(Ii)*2*RS-1;
S11i = imag(Ii)*2*RS;
S11 = 10.*log10(S11r.^2+S11i.^2);
dBVo = 20*log10(abs(Vo));
AgVo = angle(Vo)*180/pi;
uWVo = unwrap(AgVo, 179);
% uWVo = AgVo;
% toc;
% figure(2)
if logscaleEn
    semilogx(axMag, freq, dBVo, '-r', 'LineWidth', 2);
    hold(axMag, 'on');
%     semilogx(axMag, freq, S11, '-b', 'LineWidth', 2);
%     hold(axMag, 'off');
    set(axMag, 'XScale', 'log')
else
    plot(axMag, freq, dBVo, '-r', 'LineWidth', 2);
    set(axMag, 'XScale', 'linear')
end
% hold(axMag, 'on');

% w = 2.*pi.*freq;
% s = 1i.*w;
% % Hs = 1/(2*(-0.5367)).*(s.^2+1.5376)./(s.^3+2.8790.*s.^2+2.4141.*s+2.8603);
% % Hs = sqrt((16+24.*s.^2+9.*s.^4)./(16+24.*s.^2+9.*s.^4-s.^6));
% % s = s./(2.*pi);
% % Hs = (0.07)^2*((76.498+3.457e-02.*s+s.^2)./(7.4026e-01+2.0795.*s+2.93.*s.^2+2.421.*s.^3+s.^4));
% % Hs = ((76.498+3.457e-02.*s+s.^2)./(7.4026e-01+2.0795.*s+2.93.*s.^2+2.421.*s.^3+s.^4));
% Hs = 0.0002.*((53.58+0.000939.*s+s.^2)./(1.019+2.0126.*s+2.00638.*s.^2+s.^3)).^2;
% Hs_dB = 10*log10(abs(Hs));
% semilogx(axMag, freq, Hs_dB, '-m', 'LineWidth', 1);
if ~isempty(IdealFreq)
    if logscaleEn
        semilogx(axMag, IdealFreq, IdealMag, '--b', 'LineWidth', 0.1);
    else
        plot(axMag, IdealFreq, IdealMag, '--b', 'LineWidth', 0.1);
    end
end
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
hold(axPhase, 'on');
if ~isempty(IdealFreq)
    semilogx(axPhase, IdealFreq, IdealPhase, '--b', 'LineWidth', 0.1);
end
hold(axPhase, 'off');
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