%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-2-3(yyyy-mm-dd)
% 贝塞尔多项式推导
%--------------------------------------------------------------------------
function [N2, D2, ND] = fun_bessel_thomson_polynomial(m, varargin)
if nargin>1
    disp_en = varargin{1};
else
    disp_en = 0;
end
% m = 9;
N0 = zeros(1,m+1);
N1 = zeros(1,m+1);
N2 = zeros(1,m+1);
D0 = zeros(1,m+1);
D1 = zeros(1,m+1);
D2 = zeros(1,m+1);
N0(2) = vpa(0);
N1(2) = vpa(1);
D0(1) = vpa(1);
D1(1) = vpa(1);
if m == 1
    N2 = N1;
    D2 = D1;
    ND = N2+D2;
    ND = abs(ND);
    return;
end
for ii=2:m
    N2 = (2*ii-1).*N1 - [0,0,N0(1:m-1)];
    D2 = (2*ii-1).*D1 - [0,0,D0(1:m-1)];
    if disp_en
        fprintf('m=%d: ', ii);
        fprintf('N2 = ');
        st = 0;
        k  = 0;
        for jj=1:m+1
            if N2(jj) ~= 0
                if st ~= 0 && N2(jj)>0 % add +
                    fprintf(' + ');
                end
                if jj==1
                    fprintf('%d', N2(jj));
                elseif jj == 2
                    fprintf('%d*x', N2(jj));
                else
                    fprintf('%d*x^%d', N2(jj), jj-1);
                end
                st = 1;
                k = jj;
            end
        end
        fprintf('\n');
        fprintf('          ');
        for jj=1:k
            fprintf('-');
        end
        fprintf('-----\n');
        fprintf('     D2 = ');
        st = 0;
        for jj=1:m+1
            if D2(jj) ~= 0
                if st ~= 0 && D2(jj)>0 % add +
                    fprintf(' + ');
                end
                if jj==1
                    fprintf('%d', D2(jj));
                elseif jj == 2
                    fprintf('%d*x', D2(jj));
                else
                    fprintf('%d*x^%d', D2(jj), jj-1);
                end
                st = 1;
            end
        end
        fprintf('\n');
        fprintf('\n');
    end
    D0 = D1;
    D1 = D2;
    N0 = N1;
    N1 = N2;
end
ND = N2+D2;
ND = abs(ND);
end