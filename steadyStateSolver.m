function [T_final, psi_final, a_final, C_final, filename, new_b] = steadyStateSolver(r, z, t, T_n, psi_nMinusHalf, a_n, C_n, cartesian, bc_type, constants, T_coeffs, psi_coeffs, options)

show_intermediate_graphs = options(1);
steady_state_condition = options(2);

[dt, t_num] = gridProperties(t);
[~, ~, r_num, ~] = meshGridProperties(r, z);

Rm = constants('Rm');

if bc_type == 'mushyLayer'
    H = constants('H');
    R = constants('R');
    b = constants('b');
    a_fixed = constants('a_fixed');
    new_b = b; %Only use this if we need to change b because a > b
end

averaged = 0;


steadyState = false;

if show_intermediate_graphs
    clf;
    figure;
    hold on;
end

for k=1:t_num
    
    psi_bc_parameters = containers.Map({'T_n', 'C_n', 'a_n'},{T_n, C_n, a_n});
    T_bc_parameters = []; % Not used at the moment
    
    %Calculate psi_n
    f_n = poissonRHS(bc_type, r, z, T_n, constants);
    psi_n = poissonSolverTwoDimensionsAxisymm(f_n, psi_nMinusHalf, r, z, bc_type, psi_bc_parameters, constants, psi_coeffs);
    
    %Calculate T_n+1/2
    T_nPlusHalf = heatSolverADI(r, z, t, T_n, psi_n, cartesian, bc_type, T_bc_parameters, constants , T_coeffs);
    
    %Calculate psi_n+1/2
    f_nPlusHalf = poissonRHS(bc_type, r, z, T_nPlusHalf, constants);
    
    psi_bc_parameters('T_n') = T_nPlusHalf;
    psi_nPlusHalf = poissonSolverTwoDimensionsAxisymm(f_nPlusHalf, psi_n, r, z, bc_type, psi_bc_parameters, constants, psi_coeffs);
    
    %Calculate T_n+1
    T_nPlusOne = heatSolverADI(r, z, t, T_nPlusHalf, psi_nPlusHalf, cartesian, bc_type, T_bc_parameters, constants, T_coeffs);
    
    if show_intermediate_graphs
        clf;
        %mesh(r, z, psi_nPlusHalf);
        
        [u_r, u_z] = axiVelocitiesFromPsi(r, z, psi_nPlusHalf);
        streamslice(r, z, u_r.', u_z.');
        hold on;
       % axis([0 1 0 1]);
        contour(r, z, T_nPlusOne);
        xlabel('r');
        ylabel('z');
        
        hold off
        
        c = colorbar();
        c.Label.String = 'Temperature';
        
        %plot(a_n, z);
        %xlabel('a');
        %ylabel('z');
        
        
        drawnow;
       
    end    
    
    if strcmp(bc_type, 'mushyLayer')
        %C_nPlusOne = concentration(r, z, psi_nMinusHalf, T_n, C_n, a_n, constants);
        [a_nPlusOne, a_extrap] = extrapolate_a(r, z, t, a_n, constants, psi_nMinusHalf, T_n, C_n);
        
        %Keep track of how a(t) changes, and if it gets stuck at a=b then
        %stop iterating and increase b slightly
        a_values(k) = a_nPlusOne(1);
        a_extrap_values(k) = a_extrap(1);
        
        if k > 4 && ...
                a_values(k) == a_values(k-1) && ...
                a_values(k-1) == a_values(k-2) && ...
                a_values(k) == b
            
            new_b = b*1.005;
            fprintf('Increasing b to %1.5f \n', new_b);
            break; %quit iterating
        end
        
        
    end
    
    % Check for steady state
    max_T_diff(k) = max(max(abs(T_nPlusOne - T_n)))/dt;
    max_psi_diff(k) = max(max(abs(psi_nPlusHalf - psi_nMinusHalf)))/dt;
    
    if strcmp(bc_type, 'mushyLayer')
        max_a_diff(k) = max(max(abs(a_nPlusOne - a_n)))/dt;
    end

    
    %max_diff = max(max(abs((T_nPlusOne - T_n)./T_n))) + ...
    %    max(max(abs((psi_nPlusHalf - psi_nMinusHalf)./psi_nMinusHalf)));
    
    if (abs(max_T_diff(k))) < steady_state_condition &&  ...
        (abs(max_psi_diff(k))) < steady_state_condition
    
        if ~(strcmp(bc_type, 'mushyLayer') && (abs(max_a_diff(k))) > steady_state_condition)
            steadyState = true;
            break;
        end
    end
    
    % Check for oscillations
    if strcmp(bc_type, 'mushy_layer') && ...
            (averaged < 1 && k > 10) && (...
        (max_T_diff(k) > max_T_diff(k-1) && ...
        max_T_diff(k-1) < max_T_diff(k-2) && ...
        max_T_diff(k-2) > max_T_diff(k-3)) || ...
        (max_a_diff(k) > 0 && max_a_diff(k-1) < 0))
        
    
        % If we have oscillations, try averaging
        T_n = (T_nPlusOne + T_n)/2;
        psi_nMinusHalf = (psi_nMinusHalf + psi_nPlusHalf)/2;
        a_n = (a_n + a_nPlusOne)/2;
        
        averaged = averaged+1; %only average once
        fprintf('Averaging to try and remove oscillations \n');
        
    else
        % Proceed normally
        % Get ready for next loop
        T_n = T_nPlusOne;
        psi_nMinusHalf = psi_nPlusHalf;
        
        if strcmp(bc_type, 'mushyLayer')
            a_n = a_nPlusOne;
        end
    end
        
    

    if strcmp(bc_type, 'mushyLayer')
        fprintf('Completed %d out of %d time loops, max(dT/dt) = %1.10f, max(dpsi/dt) = %1.10f, max(da/dt) = %1.10f \n', k, t_num, max_T_diff(k), max_psi_diff(k), max_a_diff(k));
    else
        fprintf('Completed %d out of %d time loops, max(dT/dt) = %1.10f, max(dpsi/dt) = %1.10f \n', k, t_num, max_T_diff(k), max_psi_diff(k));
    
    end
    
 end


T_final = T_nPlusOne;

% Calculate psi for the final T state
f_nPlusOne = poissonRHS(bc_type, r, z, T_nPlusOne, constants);
psi_bc_parameters('T_n') = T_nPlusOne;
psi_final = poissonSolverTwoDimensionsAxisymm(f_nPlusOne, psi_nPlusHalf, r, z, bc_type, psi_bc_parameters, constants, psi_coeffs);

if strcmp(bc_type, 'mushyLayer')
    a_final = a_nPlusOne;
    C_final = C_n;
    
    a = a_final;    C = C_final;
end

theta = T_final;   psi = psi_final;     

if steadyState == true
    q = calculateQ(r, z, theta, psi);
    vars = filenameVars(constants, r_num, bc_type);
    filename = strcat('\data\',bc_type, 'PrevSteadyState',vars,'.mat');
    if strcmp(bc_type, 'mushyLayer')
        save([pwd filename],'r','z','theta','psi', 'a', 'C', 'constants', 'q');
    else
        save([pwd filename],'r','z','theta','psi', 'constants');
    end
    fprintf('Saved as %s\n', filename);
else
    filename = false;
    fprintf('Did not reach steady state - no file saved \n');
end

end

function f_n = poissonRHS(bc_type, r, z, T_n, constants)

[dr, dz, ~, ~] = meshGridProperties(r, z);
[dT_dr, ~] = gradient2order(T_n, dr, dz);

switch(bc_type)
    case 'mushyLayer'
        f_n = -constants('Rm') * dT_dr;
    case 'heatedWire'
        f_n = constants('Rm') * dT_dr;
    otherwise
        f_n = 0;
end

end
        