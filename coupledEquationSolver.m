function [T_nPlusOne, psi_nPlusHalf]  = coupledEquationSolver(r, z, t, T_n, psi_nMinusHalf, C_n, a_n, cartesian, bc_type, constants, Coeffs)

[dr, ~] = gridProperties(r);
[dz, ~] = gridProperties(z);

Rm = constants('Rm');

psi_bc_parameters = containers.Map({'T_n', 'C_n', 'a_n'},{T_n, C_n, a_n});
T_bc_parameters = []; % Not used at the moment

%Calculate psi_n
[~, dT_dr] = gradient(T_n, dz, dr);
f_n = Rm * dT_dr;
psi_n = poissonSolverTwoDimensionsAxisymm(f_n, psi_nMinusHalf, r, z, bc_type, psi_bc_parameters, constants);

%Calculate T_n+1/2
T_nPlusHalf = heatSolverADI(r, z, t, T_n, psi_n, cartesian, bc_type, T_bc_parameters, constants , Coeffs);

%Calculate psi_n+1/2
[~, dT_dr] = gradient(T_nPlusHalf, dz, dr);
f_nPlusHalf = Rm * dT_dr;

psi_bc_parameters('T_n') = T_nPlusHalf;
psi_nPlusHalf = poissonSolverTwoDimensionsAxisymm(f_nPlusHalf, psi_n, r, z, bc_type, psi_bc_parameters, constants);

%Calculate T_n+1
T_nPlusOne = heatSolverADI(r, z, t, T_nPlusHalf, psi_nPlusHalf, cartesian, bc_type, T_bc_parameters, constants, Coeffs);

end