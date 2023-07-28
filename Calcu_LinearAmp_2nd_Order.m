%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-13(yyyy-mm-dd)
% LinearAmplitude滤波器综合，参数计算：2阶
%--------------------------------------------------------------------------

RL = 1;
% R1 = 0.5;
% Rs = 0.3;
syms L1 C2 Z1 Z2 s
R1 = RL/2;%(sqrt(5)-1)/2*RL;
Rs = R1*RL/(R1+RL);
Z1 = RL*(1/(s*C2))/(RL+(1/(s*C2)));
Z2 = Z1+s*L1;
Z3 = Z2*R1/(R1+Z2);
Z4 = Z3+Rs;

f = Z2/Z1*Z4/Z3;
f_s = collect(f, s);
coeff_f = coeffs(f_s, s);
% coeff_target = 4/1.488*[1.488,0.976,0.488];
% RLP = RL*R1/(RL+R1);
coeff_target = fliplr([1,1.222,1.746])./1.746.*(2);
% coeff_target = [1.220, 0.853, 0.699];
eqns = arrayfun(@(a, b) a == b, coeff_f(2:end), coeff_target(2:end));
sol = solve(eqns, L1, C2);
C2 = vpa(sol.C2,4);
L1 = vpa(sol.L1,4);
fprintf('=========2nd Order Linear Amplitude Filer ==========\n');
fprintf('Rs = %s Ohm\n', Data2Enginner(Rs, '0.3'));
fprintf('R1 = %s Ohm\n', Data2Enginner(R1, '0.3'));
fprintf('RL = %s Ohm\n', Data2Enginner(RL, '0.3'));
fprintf('C2 = %s F\n', Data2Enginner(double(C2(1)), '0.3'));
fprintf('L1 = %s H\n', Data2Enginner(double(L1(1)), '0.3'));

% vpa(sol.R1,4)

% syms R1 Rs Rsp Rs0 RL

% Res = solve(R1*Rs/(R1+Rs)==Rsp, R1*(Rs0+RL)==Rs*(RL+R1)+R1*RL, R1, Rs);
% vpa(Res.R1)
% vpa(Res.Rs)
% syms R1 Rs RL Rs1 Rs0 Rp1 Rp2
% Res = solve(Rs*R1/(Rs+R1)==Rp2*(Rs1+Rs0*Rp1/(Rs0+Rp1))/(Rp2+(Rs1+Rs0*Rp1/(Rs0+Rp1))),...
%     R1/(R1+Rs)==Rp2/(Rs1+Rp2)*(Rp1*(Rs1+Rp2)/(Rp1+Rs1+Rp2))/(Rs1+(Rp1*(Rs1+Rp2)/(Rp1+Rs1+Rp2))),...
%     Rs0 == Rp1*(Rs1+Rp2*RL/(Rp2+RL))/(Rp1+(Rs1+Rp2*RL/(Rp2+RL))), Rp1, Rp2, Rs1);
% Rp1 = Res.Rp1
% Rp2 = Res.Rp2
% Rp3 = Res.Rp3

