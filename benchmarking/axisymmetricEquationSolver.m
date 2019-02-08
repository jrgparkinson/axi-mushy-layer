%{
Specify initial conditions and grid and then use the SteadyStateSolver
to solve the heat and momentum equations
%}

function axisymmetricEquationSolver

clear; close all;

show_intermediate_graphs = false;
show_final_graphs = false;

% Specify the boundary conditions for this problem. 
bcs = 'heatedWire';

%Options for heated wire analytic solution
reltol = 1e-7;
abstol = 1e-8;
bvp_options = bvpset('RelTol',reltol, 'AbsTol', abstol);

steady_state_condition = 1e-6;

a_nc = 0.2;
Rm = 10; %Rayleigh number.
T_inf = 1;
alpha_m = 1;

%Coeffs contains coefficients of terms:
%r T_t - A*psi_z*T_r + B*psi_r*T_z = C*r*T_rr + D*r*T_zz  + E*(r*T_r)_r
%Coeffs = [A, B, C, D, E];
Coeffs = [1, 1, 0, alpha_m, alpha_m];

r_min = sqrt(a_nc/Rm);
r_max = r_min +1;
z_min = 0;
z_max = z_min + 0.5;

% loop for different grid spacings
grid_points = [15 20 25 30];
drs = (r_max-r_min) ./ (grid_points-1);

% Run for large number of timesteps to ensure we hit steady state
t_num = 1e4;
dt = 1e-3;
t = 0:dt:(t_num*dt);

tic;

for grid_i = 1:numel(grid_points)
    r_num = grid_points(grid_i);
    z_num = r_num;
    
    dr = (r_max-r_min)/(r_num-1);
    dz = (z_max-z_min)/(z_num-1);
    
    r=r_min:dr:r_max;   z=z_min:dz:z_max;
    
    % Get analytic solution
    [T_analytic, psi_analytic, ~, ~] = heatedWireAnalyticSolution(r, z, a_nc, T_inf, Rm, bvp_options, alpha_m);
    
    % Set initial conditions
    T_n = T_analytic;
    psi_nMinusHalf = psi_analytic;
    
    % Solve for steady state
    options = [show_intermediate_graphs, steady_state_condition];
    [T_final, psi_final] = steadyStateSolver(r, z, t, T_n, psi_nMinusHalf, false, bcs, Rm, Coeffs, options);
    
    % Data analysis
    frac_error = abs((T_final - T_analytic) ./ (T_analytic + 1));
    psi_frac_error = abs( (psi_final - psi_analytic) ./ (psi_analytic + 1));
    
    if show_final_graphs
        %mesh(T_n);
        %title('initial');% Plot the initial solution.
        %pause;
        var_string = strcat('steady state, r T_t = \alpha_m (r T_r)_r + \psi_x T_r - \psi_r T_x, \Delta t = ', num2str(dt));
        mesh(z, r, psi_final);
        xlabel('z');
        ylabel('r');
        zlabel('\psi');
        title(strcat('Numeric solution, ', var_string));% Plot the computed solution.
        pause;
        mesh(z, r, psi_analytic);
        title(strcat('Analytic solution, ', var_string));
        xlabel('z');
        ylabel('r');
        zlabel('\psi');
        pause;
        
        mesh(z, r, psi_frac_error)          % Mesh plot of the error
        xlabel('z');
        ylabel('r');
        zlabel('(\Delta \psi)/\psi_{analytic}');
        title(strcat('Fractional error, ', var_string));
        pause;
    end
    
    av_frac_errors(grid_i) = nanmean(nanmean(frac_error));
    max_frac_errors(grid_i) = max(max(abs(frac_error.*isfinite(frac_error))));
end

timeSpent = toc;
fprintf('\nTime spent: %2.2f \n', timeSpent);

% If there's something worth plotting, plot it
if numel(drs) > 1
    %var_string = strcat('1 timestep, T_t = T_{zz} + (r*T_r)_r/r + psi_z*T_r/r - psi_r*T_z/r, \Delta t = ', num2str(dt));
    var_string = strcat('steady state, r T_t = \alpha_m (T_{zz} + (r T_r)_r) + \psi_x T_r - \psi_r T_x, \Delta t = ', num2str(dt));
    
    drs2 = drs.*drs;
    
    plot(drs2, av_frac_errors, '--*');
    xlabel('(\Delta r)^2');
    ylabel('Fractional error');
    title(strcat('Average fractional error, ', var_string));
    axis([0 max(drs2)*1.2 0 max(av_frac_errors)*1.2]);
    str = ['RelTol = ', num2str(reltol), ', AbsTol = ', num2str(abstol)];
    text(0,max(av_frac_errors),str);
    pause;
    
    plot(drs2, max_frac_errors, '--*');
    xlabel('(\Delta r)^2');
    ylabel('Fractional error');
    title(strcat('Maximum fractional error, ', var_string));
    
    fprintf('dt: %1.10f, steady state: %1.10f \n', dt, steady_state_condition);
    disp(drs.*drs);
    disp(av_frac_errors);
    
end

end



