%{
Solves one timestep of
T_t = T_rr + T_zz + (r*T_r)_r/r + psi_z*T_r/r - psi_r*T_z/r
    = T_rr + T_zz + (r*T_r)_r/r - u_r*T_r/r - u_z*T_z/r
where u = (-psi_z, psi_r)
depending on which coefficients are passed to the function
several opions:
- cartesian removes 1/r factors
- bcs specifies which boundary conditions we use. See function at bottom of file 
 to see what the options are.
%}

function T_nPlusOne = heatSolverADI(r, z, t, T_n, psi_n, cartesian, bcs, bc_parameters, constants, coeffs)

[dr, dz, r_num, z_num] = meshGridProperties(r, z);
[dt, ~] = gridProperties(t);

Coeff_T_r = coeffs('T_r');          Coeff_T_z= coeffs('T_z');
Coeff_T_z_frame = coeffs('T_z_frame');
Coeff_T_rr= coeffs('T_rr');         Coeff_T_zz= coeffs('T_zz');
Coeff_T_rf_r= coeffs('T_rf_r');

%Discretization values
T_zz = Coeff_T_zz * 1/(dz*dz);
T_rr = Coeff_T_rr * 1/(dr*dr);
T_r = Coeff_T_r * 1/(2*dr);
T_z = Coeff_T_z * 1/(2*dz);
T_z_frame = Coeff_T_z_frame * 1/(2*dz);
T_rf_r = Coeff_T_rf_r * 1/(dr*dr);

[psi_r, psi_z] = gradient(psi_n, dr, dz);

%--- sweep in r-direction --------------------------------------

T_nPlusHalf = getBoundaryConditions(bcs, r, z, T_n, psi_n, psi_z, psi_r, bc_parameters, constants);

for j = 2:z_num-1,                             
     
    A = eye(r_num);    rhs = T_nPlusHalf(j, :);
    
    for i=2:r_num-1,
        rmh = r(1, i) - dr/2;
        rph = r(1, i) + dr/2;
        r_i = r(1, i);
        if cartesian
            r_i = 1;    rmh=1;  rph=1;
        end
        
        p_z = psi_z(j, i);   p_r = psi_r(j, i);
        
        %My attempt at making A
        A(i, i-1) = -T_rr + p_z*T_r/r_i - T_rf_r*rmh/r_i;
        A(i, i) = 2/dt + 2*T_rr + T_rf_r*2*r_i/r_i;
        A(i, i+1) = -T_rr - p_z*T_r/r_i - T_rf_r*rph/r_i;
        
        rhs(i) = (T_zz - (-p_r)*T_z/r_i - T_z_frame)  *T_n(j-1, i) ...
            + (2/dt -2*T_zz)       *T_n(j, i) ...
            + (T_zz + (-p_r)*T_z/r_i + T_z_frame)       *T_n(j+1, i);
    end
    
    T_nPlusHalf(j, :) = A\(rhs.');
    
end                                 

%-------------- loop in z -direction --------------------------------

T_nPlusOne = getBoundaryConditions(bcs, r, z, T_n, psi_n, psi_z, psi_r, bc_parameters, constants);

for i = 2:r_num-1,
    
    rmh = r(1, i) - dr/2;
    rph = r(1, i) + dr/2;
    r_i = r(1, i);
    if cartesian
        r_i = 1;    rmh=1;  rph=1;   %Turns T_flux into T_rr
    end
    
    A = eye(z_num);     rhs = T_nPlusOne(:, i).';
    
    for j=2:z_num-1,
        
        p_r = psi_r(j, i);       p_z = psi_z(j, i);
        
        A(j,j-1) = -T_zz + (-p_r)*T_z/r_i + T_z_frame;
        A(j,j) = 2/dt + 2*T_zz;
        A(j,j+1) = -T_zz - (-p_r)*T_z/r_i - T_z_frame;
        
        rhs(j) = (T_rr - p_z*T_r/r_i + T_rf_r*rmh/r_i) * T_nPlusHalf(j, i-1) ...
            + (2/dt -2*T_rr - 2*T_rf_r*r_i/r_i) *T_nPlusHalf(j, i) ...
            + (T_rr + p_z*T_r/r_i + T_rf_r*rph/r_i) * T_nPlusHalf(j, i+1);
    end
    
    T_nPlusOne(:, i) = A\(rhs.');
    
    
end                             


end

function T_nPlusHalf = getBoundaryConditions(bcs, r, z, T_n, psi_n, psi_z, psi_r, bc_parameters, constants)

switch bcs
    case 'heatedWire'
        T_nPlusHalf = heatedWireBoundaryConditions(T_n);
    case 'HRL'
        T_nPlusHalf = HRLBoundaryConditions(T_n);
    case 'mushyLayer'
        T_nPlusHalf = mushyLayerTemperatureBoundaryConditions(T_n, psi_n, r, z, constants('b'));
    otherwise
        error('heatSolverADI:bcs', 'Unknown boundary conditions specified');
end

end