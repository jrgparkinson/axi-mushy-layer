function T_nPlusHalf = heatedWireBoundaryConditions(T_n)

% Get the right grid size
T_nPlusHalf = 0*T_n;

% Dirichlet boundary conditions
T_nPlusHalf(:,1) = T_n(:, 1);
T_nPlusHalf(:,end) = T_n(:,end);
T_nPlusHalf(1,:) = T_n(1, :);
T_nPlusHalf(end,:) = T_n(end,:);

end

