 function psi_nPlusHalf = mushyLayerPsiBoundaryConditions(r, z, psi_n, bc_parameters, constants)

% Create empty grid of the correct size
psi_nPlusHalf = psi_n;

T_n = bc_parameters('T_n');
C_n = bc_parameters('C_n');
a_n = bc_parameters('a_n');

b = constants('b');
Rm = constants('Rm'); 
Da = constants('Da');

[dr, dz, r_num, z_num] = meshGridProperties(r, z);


%Left boundary, r=b, complicated - chimney
%See latex notes for details of the numerical approximation
A = a_n(:, 1).^4/(16 * Da);
denom = 1 + (1.5/dr) * (0.5*b + A/b);

% B is a function of z only
% 0th order approx:
%B = A*Rm.*(T_n(:, 1) - C_n(:, 1)); 
% 1st order approx:
B = A*Rm.*(T_n(:, 1) - C_n(:, 1)) + (1/6) * (b^3 - a_n(:, 1).^3) * Rm .* (-3*T_n(:, 1) + 4*T_n(:, 2) - T_n(:, 3))/(2*dr);

psi_nPlusHalf(:, 1) = (1./denom) .* ((0.5*b + A/b).* ( (4 * psi_n(:, 2) - psi_n(:, 3))/(2*dr) ) + B   );


%Top boundary, z=H, psi_z = 0
psi_nPlusHalf(end, :) = (1/3) * (4 * psi_n(end-1, :) - psi_n(end-2, :));

%Bottom boundary, z=0, psi = 0
psi_nPlusHalf(1, :) = 0;

%Right boundary, r=H, psi = 0
psi_nPlusHalf(:, end) = 0;

end