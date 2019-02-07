 %{
Solve for the steady state of a mushy layer in an axisymmetric geometry
Boundary conditions specified in

- mushyLayerTemperatureBoundaryConditions
- poisson.

Initial conditions from: mushyLayerInitialState
%}

function axisymmetricMushyLayer

clear; close all;

% Decide whether to plot various graphs
show_intermediate_graphs = true;

steady_state_condition = 10^-5;

options = [show_intermediate_graphs, steady_state_condition];

% Specify the parameters for the model
Rms = 60; 
Da = 5*10^-5;
Hs = [0.25, 0.26, 0.27];
Rs = 0.25; %:0.006:0.352;
a_fixed = -0.0325; %Set negative to use relaxation
relaxation_lambda = 0.002; %usually 0.002
b = Rs/10; % Initial guess at chimney width
chimney_method = 3;  %options: 1 = f_extrap, 2 = Q_contour, 3 = cfit
r_num = 60;  z_num = 60;

file_to_load = false;

% Specify the boundary conditions for this problem
bc_type = 'mushyLayer';
cartesian = false;

% Run for large number of timesteps to ensure we hit steady state
t_num = 2000;
dt = 1e-1; %usually 1e-1
t = 0:dt:((t_num-1)*dt);

%Determine which equations are being solved
T_coeffs = containers.Map({'T_r', 'T_z', 'T_z_frame', 'T_rr', 'T_zz', 'T_rf_r'}, [1, 1, 1, 0, 1, 1]);
psi_coeffs = containers.Map({'psi_r', 'psi_zz'}, [1, 1]);

for Rm_i = 1:numel(Rms)
    Rm = Rms(Rm_i);

    %Iterate over different R
    for R_i = 1:numel(Rs)
        R = Rs(R_i);
        fprintf('Starting a new run with R = %1.4f \n', R);

        %Iterate over different H
        for H_i = 1:numel(Hs)
            H = Hs(H_i);
            fprintf('Starting a new run with H = %1.4f \n', H);

            %{
            if R_i > 1 && R < Rs(R_i - 1)
                %If we're decreasing R, decrease b to the last a
                b = a_final(1);
                fprintf('Have decreased b to %1.6f\n', b);
            elseif R_i > 1
                b = a_final(1);
            %}
            if H_i > 1 || Rm_i > 1 || R_i > 1
                b = a_final(1); %Trying to ensure b is as close to the new a as possible
            end

            % Set initial conditions
            % See if we already have data for these parameters
            flst=dir([pwd '/data/*.mat']);
            filelist={flst.name};
            
            existing_data = find_file_with_vars(R,H,Rm, r_num, filelist);
            if existing_data
                file_to_load = existing_data;

                % Extract b from this data
                b_format = 'b(\d\.\d+)';
                [matches] = regexp(file_to_load, b_format, 'match');
                b = str2num(strrep(matches{1}, 'b', ''));
            else
                file_to_load = false;
            end

            %file_to_load = '/data/mushyLayerPrevSteadyState_relaxed_Points40Rm60H0.25R0.232b0.032101.mat';
            %force this for now
            %b=0.032;

            %Create the grid
            r_min = b;    r_max = R;
            z_min = 0;    z_max = H;

            dr = (r_max-r_min)/(r_num-1);   dz = (z_max-z_min)/(z_num-1);

            [r, z] = meshgrid(r_min:dr:r_max, z_min:dz:z_max);

            % Put all the constants in one object to make it easier to pass them around
            constants = containers.Map({'Rm','Da', 'R', 'H', 'relaxation_lambda', 'chimney_method', 'b', 'a_fixed'}, ...
                [Rm, Da, R, H, relaxation_lambda, chimney_method, b, a_fixed]);

            % Get the initial data on this grid
            [T_initial, psi_initial, a_initial, C_initial] = mushyLayerInitialState(r, z, constants, file_to_load);

            % Solve for steady state
            [T_final, psi_final, a_final, C_final, file_to_load, new_b] = steadyStateSolver(r, z, t, T_initial, psi_initial, a_initial, C_initial, cartesian, bc_type, constants, T_coeffs, psi_coeffs, options);

            % Check to see if we had an issue with b being < a, if so try again
            while new_b ~= b
                b = new_b;  constants('b') = b;

                [T_final, psi_final, a_final, C_final, file_to_load, new_b] = steadyStateSolver(r, z, t, T_final, psi_final, a_final, C_final, cartesian, bc_type, constants, T_coeffs, psi_coeffs, options);
            end

            % If we couldn't get a steady state, move to the next R value and
            % vary H along there
            if file_to_load == false
               break; 
            end

            % Clear up ready for the next run?
            
            theta_inf = calculateThetaInfinity(r,z,T_final, psi_final);
            if theta_inf < 1.4
                H_i = numel(Hs);
            end
            
        end

    end
end

fprintf('Finished!\n');
end
