%Solve the HRL problem

function HRL

clear; close all;

bc_type = 'HRL';

% Decide whether to plot various graphs
show_intermediate_graphs = false;
steady_state_condition = 1*10^-5;
options = [show_intermediate_graphs, steady_state_condition];

%Parameters
%Rms = [10 4*pi^2 50 100 200 250 300];
%x_nums = [24 24 32 32 48 48 48];
Rms = [10 4*pi^2 50 100];
x_nums = [24 24 32 32];
A = 1; %aspect ratio

Nus = nan*ones(numel(Rms));

L = 1; H = A*L;
z_min = 0;  z_max = H;
x_min = 0;  x_max = L; 

t_num = 2000;
dt = 1e-2;
t = 0:dt:((t_num-1)*dt);

for i = 1:numel(Rms)
    Rm = Rms(i);
    x_num = x_nums(i); z_num = x_num;
    
    dx = L/(x_num - 1);   dz = (z_max-z_min)/(z_num-1);
    
    [x, z] = meshgrid(x_min:dx:x_max, z_min:dz:z_max);
    
    constants = containers.Map({'Rm','Da', 'R', 'H', 'relaxation_lambda', 'chimney_method', 'b', 'a_fixed'}, ...
        [Rm, 0, 0, 0, 0, 0, 0, 0]);
    
    %Determine which equations are being solved
    T_coeffs = containers.Map({'T_r', 'T_z', 'T_z_frame', 'T_rr', 'T_zz', 'T_rf_r'}, [1, 1, 0, 1/Rm, 1/Rm, 0]);
    
    
    %Try and find a file for this data
    filename = strcat('\data\',bc_type,'\', 'steadyStateRm', num2str(Rm),'x_num',num2str(x_num),'.mat');
    if(exist([pwd filename], 'file') == 2)
        load([pwd filename]);
    else
        % Get the initial data on this grid
        [T_initial, psi_initial] = HRLInitialState(x, z, constants);

        % Solve for steady state
        [T_final, ~] = steadyStateSolverCartesian(x, z, t, T_initial, psi_initial, T_coeffs, constants, options, bc_type);
    end
    
    mesh(x, z, psi_final); pause;
    mesh(x, z, T_final); pause;
    
    % Calculate the nusselt number
    [~, T_z] = gradient2order(T_final, dx, dz);
    Nu = -trapz(x(1, :),T_z(1, :));
    Nus(i) = Nu;
    
    
    
end

plot(Rms, Nus);
