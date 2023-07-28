%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-13(yyyy-mm-dd)
% LinearAmplitude滤波器综合，参数计算：4阶
%--------------------------------------------------------------------------
n = 2;
Ap = 10*log10(2);% 3dB 衰减
RL = 1;
% R1 = 0.5;
% Rs = 0.3;
% 从负载往原依次为RL, C1, L2, C3, L4, ..., R1, Rs
syms s C1 L2 C3 L4 C5 L6
R1 = (sqrt(5)-1)/2*RL;
Rs = R1*RL/(R1+RL);
Z1 = RL*(1/(s*C1))/(RL+(1/(s*C1)));
Z2 = Z1+s*L2;
% Z3 = Z2*(1/(s*C3))/(Z2+(1/(s*C3)));
% Z4 = Z3+s*L4;
Zs1 = Z2*R1/(R1+Z2);
Zs = Zs1+Rs;

f = Z2/Z1*Zs/Zs1;

IL   = 10^(-Ap/20);
[P, w1, absND] = fun_linear_amp_polynomial(n, IL);
ND     = (funRecursionPoly(n, P));
ND2 = eval(ND);
ND3 = real(ND2);
% absND = abs(eval(absND));
f_s = collect(f, s);
coeff_f = coeffs(f_s, s);
% coeff_target = 4/1.488*[1.488,0.976,0.488];
% RLP = RL*R1/(RL+R1);
coeff_target = (ND3)./ND3(1).*2;% 2倍增益
% coeff_target = [1.220, 0.853, 0.699];
eqns = arrayfun(@(a, b) a == b, coeff_f(2:end), coeff_target(2:end));
sol = solve(eqns)
C1 = vpa(sol.C1,4);
L2 = vpa(sol.L2,4);
% C3 = vpa(sol.C3,4);
% L4 = vpa(sol.L4,4);
fprintf('=========2nd Order Linear Amplitude Filer ==========\n');
fprintf('Rs = %s Ohm\n', Data2Enginner(Rs, '0.3'));
fprintf('R1 = %s Ohm\n', Data2Enginner(R1, '0.3'));
fprintf('RL = %s Ohm\n', Data2Enginner(RL, '0.3'));
fprintf('C1 = %s F\n', Data2Enginner(double(C1(1)), '0.3'));
fprintf('L2 = %s H\n', Data2Enginner(double(L2(1)), '0.3'));
% fprintf('C3 = %s F\n', Data2Enginner(double(C3(1)), '0.3'));
% fprintf('L4 = %s H\n', Data2Enginner(double(L4(1)), '0.3'));

% vpa(sol.R1,4)