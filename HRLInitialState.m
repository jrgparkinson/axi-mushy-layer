function [T_initial, psi_initial] = HRLInitialState(x, z, constants)

Rm = constants('Rm');

[~, ~, x_num, z_num] = meshGridProperties(x,z);

epsilon = 0.1;

T = zeros(x_num, z_num);            psi = zeros(x_num, z_num);
T_analytic = zeros(x_num, z_num);   psi_analytic = zeros(x_num, z_num);

T = 1-z + epsilon*cos(pi*x).*sin(pi*z);
delta = Rm * epsilon / (2*pi);
psi = delta * sin(pi*x) .* sin(pi*z);

%sigma = 2*pi^2 - Rm/2;
%T_analytic(i, j) = 1 - z + epsilon*cos(pi*x)*sin(pi*z)*exp(sigma*t_max);
%psi_analytic(i, j) = delta * sin(pi*x) * sin(pi*z)*exp(sigma*t_max/2);

%Enforce B.C's
T(:, z_num) = 0;
T(:, 1) = 1;


T_initial = T;
psi_initial = psi;

end