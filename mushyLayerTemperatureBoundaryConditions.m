function T_nPlusHalf = mushyLayerTemperatureBoundaryConditions(T_n, psi_n, r, z, b)

[dr, dz, r_num, z_num] = meshGridProperties(r, z);

% Get the right grid size
T_nPlusHalf = T_n;


%Chimney boundary condition, r=b, psi_q theta_z(b) = b theta_r(b)
%See latex notes for details of numerical approximation

psi_q = psi_n(:, 1) - 0.5 * b^2;
for j = 2:(z_num-1)
    T_nPlusHalf(j, 1) = (1/3) * (  -(psi_q(j) * dr)/(b * dz)  * ( T_n(j+1, 1) - T_n(j-1, 1) ) + 4 * T_n(j, 2) - T_n(j, 3) );
end

%Treat first and last points differently in order to maintain second order
%derivatives
T_nPlusHalf(1, 1) = (1/(3 * (1 - psi_q(1) * dr/(b*dz) ))) * (  (psi_q(1) * dr)/(b * dz)  * (-4*T_n(2, 1) + T_n(3, 1)) + 4 * T_n(1, 2) - T_n(1, 3));
T_nPlusHalf(end, 1) = (1/(3 * (1 + psi_q(end) * dr/(b*dz) ))) * (  (psi_q(end) * dr)/(b * dz)  * (4*T_n(end-1, 1) - T_n(end-2, 1)) + 4 * T_n(end, 2) - T_n(end, 3) );

% Neumann boundary condition for right hand side (r=R), theta_r = 0
T_nPlusHalf(:, end) = (1/3) * (4 * T_n(:, end-1) - T_n(:, end-2));

% Dirichlet boundary conditions for top and bottom (z=0, H)
T_nPlusHalf(1, :) = -1;
T_nPlusHalf(end, :) = 0;

end