function [T_final, psi_final] = steadyStateSolverCartesian(x, z, t, T_initial, psi_initial, T_coeffs, constants, options, bc_type)

show_intermediate_graphs = options(1);
steady_state_condition = options(2);

Rm = constants('Rm');

[dt, t_num] = gridProperties(t);
[dx, dz, x_num, ~] = meshGridProperties(x, z);

T_n = T_initial;
psi_nMinusHalf = psi_initial;

cartesian = true;
steadyState = false;

for k=1:t_num
    
    T_bc_parameters = []; % Not used at the moment
    
    %Calculate psi_n
    f_n = poissonRHS(dx,dz,T_n, Rm);
    psi_n = poissonSolverTwoDimensions(f_n, psi_nMinusHalf, x, z);
    
    %Calculate T_n+1/2
    T_nPlusHalf = heatSolverADI(x, z, t, T_n, psi_n, cartesian, bc_type, T_bc_parameters, constants , T_coeffs);
    
    %Calculate psi_n+1/2
    f_nPlusHalf = poissonRHS(dx,dz,T_nPlusHalf, Rm);
    
    psi_nPlusHalf = poissonSolverTwoDimensions(f_nPlusHalf, psi_n, x, z);
    
    %Calculate T_n+1
    T_nPlusOne = heatSolverADI(x, z, t, T_nPlusHalf, psi_nPlusHalf, cartesian, bc_type, T_bc_parameters, constants, T_coeffs);

    % Check for steady state
    max_T_diff(k) = max(max(T_nPlusOne - T_n))/dt;
    max_psi_diff(k) = max(max(psi_nPlusHalf - psi_nMinusHalf))/dt;
   
    
    if (abs(max_T_diff(k))) < steady_state_condition &&  ...
        (abs(max_psi_diff(k))) < steady_state_condition
        
        steadyState = true;
        break;
    end
    
    T_n = T_nPlusOne;
    psi_nMinusHalf = psi_nPlusHalf;
  
    fprintf('Completed %d out of %d time loops, max(dT/dt) = %1.10f, max(dpsi/dt) = %1.10f \n', k, t_num, max_T_diff(k), max_psi_diff(k));
end




% Calculate psi for the final T state
f_nPlusOne = poissonRHS(dx,dz,T_nPlusOne, Rm);
psi_final = poissonSolverTwoDimensions(f_nPlusOne, psi_nPlusHalf, x, z);

T_final = T_nPlusOne;

if steadyState
    filename = strcat('\data\',bc_type,'\', 'steadyStateRm', num2str(Rm),'x_num',num2str(x_num),'.mat');
    save([pwd filename], 'T_final','psi_final','x','z','constants');
end

end


function f = poissonRHS(dx,dz,T, Rm)
    [T_x, ~] = gradient2order(T, dx, dz);
    f = - T_x;
end

    