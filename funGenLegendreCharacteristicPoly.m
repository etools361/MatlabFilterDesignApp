%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-06-11(yyyy-mm-dd)
% 勒让德滤波器系数设计
% 注意返回值低次在前，高次在后
%--------------------------------------------------------------------------
function [Ln] = funGenLegendreCharacteristicPoly(FilterOrder)
% Initialize coefficient vector
coefficients = ((1:FilterOrder).*2 - 1);

% Define coefficient for normalization
if mod(FilterOrder, 2) == 1
    % Odd polynomial order
    half_order = (FilterOrder - 1) / 2;
    norm_coefficient = 1 / (sqrt(2 * (half_order + 1) * (half_order + 1)));
else
    % Even polynomial order
    half_order = FilterOrder / 2 - 1;
    norm_coefficient = 1 / sqrt((half_order + 1) * (half_order + 2));
    % Check if the half_order is odd or even
    if mod(half_order, 2) == 1
        % Odd
        coefficients(1:2:end) = 0;
    else
        % Even
        coefficients(2:2:end) = 0;
    end
end
% Apply normalization to coefficients
coefficients = coefficients * norm_coefficient;

% Initial value for the Phi
Phi = 0;

% Generate polynomial and add to Phi
for i = 1:half_order + 1
    current_coefficient = coefficients(i);
    legendrePoly = JacobiPoly(i - 1, 0, 0); % Gen legendre polynomial
    Phi  = polyadd(Phi, current_coefficient * legendrePoly);
end
% Calculate Phi^2
Phi_squared = conv(Phi, Phi); 

% For even polynomial order, multiply Phi_squared by (x+1)
if mod(FilterOrder, 2) == 0
    Phi_squared = conv(Phi_squared, [1,1]);
end
% Calculate the integral of Phi_squared
integral_of_Phi_squared = polyint(Phi_squared, 0);

% Define the range for polyval function
lower_range = -1; 
upper_range = [2, -1];

% Calculate polyval for lower and upper ranges
D = polyval(integral_of_Phi_squared, lower_range);
U = substitute_poly(integral_of_Phi_squared, upper_range);

% Final result
Ln0 = polyadd(U, -D);
nLn = length(Ln0);
Ln = zeros(1, nLn*2-1);
Ln(3:2:end) = fliplr(Ln0(1:end-1));
Ln(1)     = Ln0(end);

end



