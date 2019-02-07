function [T_initial, psi_initial, a_initial, C_initial] = mushyLayerInitialState(r, z, constants, file_to_load)

H = constants('H');
a = constants('a_fixed');

[~, ~, r_num, z_num] = meshGridProperties(r, z);

empty_grid = 0*meshgrid(1:r_num, 1:z_num);
theta = empty_grid; psi = empty_grid;

r_interp = r;       z_interp = z;

if file_to_load
    load([pwd file_to_load])
else
    
    % Try and find some existing data
    vars = filenameVars(constants, r_num);
    filename = strcat('\data\mushyLayerPrevSteadyState',vars,'.mat');
    
    if exist([pwd filename], 'file')
        load([pwd filename]);
    else
        %Default option
        %load([pwd '\data\mushyLayerPrevSteadyStatePoints40Rm60H0.5R0.5b0.15.mat']);
        load([pwd '\data\mushyLayerPrevSteadyStatePoints40Rm60H0.25R0.25a0.0325b0.0325.mat']);
    end
end


T_initial = interp2(r, z, theta, r_interp, z_interp, 'cubic', 0);
psi_initial = interp2(r, z, psi, r_interp, z_interp, 'cubic', 0);

if exist('a') == 1 && a(1) > 0
    a_initial = interp1(z(:, 1), a, z_interp(:, 1));
else
    a_initial =  0.032232 * ones(z_num, 1);
end

C_initial = -1 + z_interp/(2*H);

end
