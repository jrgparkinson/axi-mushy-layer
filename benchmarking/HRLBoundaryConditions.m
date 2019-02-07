function T_nPlusHalf = HRLBoundaryConditions(T_n)

% Get the right grid size
T_nPlusHalf = 0*T_n;

% Dirichlet boundary conditions for z=0,1
T_nPlusHalf(:,1) = T_n(:, 1);
T_nPlusHalf(:,end) = T_n(:,end);

% Neumann boundary conditions for r=0,1
T_nPlusHalf(1, :) = (1/3) * (4 * T_n(2, :) - T_n(3, :));
T_nPlusHalf(end, :) = (1/3) * (4 * T_n(end-1, :) - T_n(end-2, :));

end

