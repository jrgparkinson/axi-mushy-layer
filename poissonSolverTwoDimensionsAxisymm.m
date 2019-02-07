%==========================================================================
% Poisson Solver
% Solves psi_rr / r + psi_zz / r - psi_r/(r^2) = f(r,z) on the r-z grid specified.
% Uses SOR, with the parameter alpha.
% Dirichlet boundary conditions given by psi_init
%==========================================================================
function psi_now = poissonSolverTwoDimensionsAxisymm(f, psi_init, r, z, bc_type, bc_parameters, constants, psi_coeffs)

%Determine the convergence requirement based on how large the values of psi
%are
convergence_criteria = 10^-8 * (max(psi_init) + 1); %add 1 incase psi_init is entirely 0
alpha = 1.8; %SOR parameter

[dr, dz, r_num, z_num] = meshGridProperties(r, z);

coeff_r = psi_coeffs('psi_r');
coeff_zz = psi_coeffs('psi_zz');

%Initializing previous iterations
psi_prev=psi_init;

%Boundary conditions
switch bc_type
    case 'mushyLayer'
        psi_now = mushyLayerPsiBoundaryConditions(r, z, psi_init, bc_parameters, constants);
    otherwise
        psi_now = psi_init;
end

%Giving initial difference between psi_now and psi_prev to start the
%iterations
psi_prev(round(r_num/2),round(z_num/2))=psi_prev(round(r_num/2),round(z_num/2))+0.01;

%Iteration loop
while(max(max(abs(psi_now-psi_prev)))>convergence_criteria)%Run this until convergence
   
    psi_prev=psi_now;

    for i=2:(r_num-1)
        r_i = r(1, i);
        rmh = r_i - dr/2;
        rph = r_i + dr/2;
        
        for j=2:(z_num-1)
            
            A = coeff_r  * 1/(rmh*dr^2);
            B = coeff_r  * 1/(rph*dr^2);
            C = coeff_zz * 1/(r_i*dz^2);
            
            denom = A + B + 2*C;
            
            psi_now(j, i)=(1-alpha)*psi_now(j, i)+ (alpha/denom) * ( - f(j, i) + ...
                A * psi_now(j, i-1) + B * psi_now(j, i+1) +  ...
                C * (psi_now(j+1, i) + psi_now(j-1, i))   );
        
        end
        
    end

end

end