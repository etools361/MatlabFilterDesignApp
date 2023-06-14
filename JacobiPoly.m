%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-05-31(yyyy-mm-dd)
% 雅克比多项式计算
% ref:Lowpass Filters Approximation Based on the Jacobi Polynomials, Nikola
% Stojanovic, Negovan Stamenkovic.
%--------------------------------------------------------------------------
function P=JacobiPoly(n,a,b)
% Coefficients P of the Jacobi polynomial
% They are stored in decending order of powers
if nargin == 1
    a=0; b=0;
elseif nargin == 2
    b=a;
end
p0 = 1;
p1 = [(a+b)/2+1,(a-b)/2];
if n == 0
    P=p0;
elseif n == 1
    P=p1;
else
    for k=2:n
        d=2*k*(k+a+b)*(2*k-2+a+b);
        A=(2*k+a+b-1)*(2*k+a+b-2)*(2*k+a+b)/d;
        B=(2*k+a+b-1)*(a^2-b^2)/d;
    %     Lowpass filters approximatin based on the jacobi polynomials 5
        C=2*(k-1+a)*(k-1+b)*(2*k+a+b)/d;
        P=conv([A B],p1)-C*[0,0,p0];
        p0 = p1;
        p1 = P;
    end
end
end