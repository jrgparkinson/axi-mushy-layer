function nusseltNumber
clear all;
close all;
show_intermediate_graph = false;

steady_state_criteria = 10^-6;

%Constants - don't change these!
Ra_crit = 4*(pi^2); %roughly 39.5
Rms = [300];
epsilon = 1e-1; %1e-5 to remove non-linear effects

SOR_parameter = 1.8; %Relaxation parameter
%Setup the grid
z_min = 0;  z_max = 1;          
x_min = 0;  x_max = 1;     L = x_max - x_min;                          
t_min = 0; 

x_nums = [48];%20 22 24 26 28 30 40 50 80
dxs = L ./ (x_nums - 1);
    
%t_num = round(x_num * u_max * (t_max - t_min) / (x_max - x_min)); %cfl condition
%dt = (t_max - t_min)/t_num;
t_max = 50;
dts = 1e-3*[10]; %usualy 10e-3 i.e. 1e-2

neumann = true;
cartesian = true;

%loop for dt
for dt_i = 1:numel(dts)
    dt = dts(dt_i);
    t_num = round(t_max/dt);
    
    mod_u_2 = 0;
    
    %Loop for dx
    for x_num_i = 1:numel(x_nums)
        x_num = x_nums(x_num_i);
        z_num = x_num;  
        x_points = linspace(x_min, x_max, x_num);
        z_points = linspace(z_min, z_max, z_num);
        
        dx = (x_max - x_min)/(x_num-1); dz = (z_max - z_min)/(z_num-1);
        
        %Loop for Rm
        for rm_i = 1:numel(Rms);
            Rm = Rms(rm_i); %Rayleigh number  
            
            T_xx = 1/(dx^2);
            T_zz = 1/(dz^2);
            T_x = 1/(2*dx);
            T_z = 1/(2*dz);
            T_xf_x = 0;

            % Create initial grids
            T = zeros(t_num, x_num, z_num);     T_analytic = zeros(x_num, z_num);
            psi = zeros(t_num, x_num, z_num);   psi_analytic = zeros(x_num, z_num);
            for i = 1:x_num
                x = x_points(i);
                for j = 1:z_num
                    z = z_points(j);
                    T(1, i, j) = 1-z + epsilon*cos(pi*x)*sin(pi*z);
                    delta = Rm * epsilon / (2*pi);
                    psi(1, i, j) = delta * sin(pi*x) * sin(pi*z);

                    sigma = 2*pi^2 - Rm/2;

                    T_analytic(i, j) = 1 - z + epsilon*cos(pi*x)*sin(pi*z)*exp(sigma*t_max);
                    psi_analytic(i, j) = delta * sin(pi*x) * sin(pi*z)*exp(sigma*t_max/2);
                end
                
                %Enforce B.C's
                T(1, i, z_num) = 0;
                T(1, i, 1) = 1;
            end

            diff = 10000; %Difference between time steps, set as something large to start

            if show_intermediate_graph
                figure;
                hold on;
            end

            %Iterate over time
            for t_n = 1:t_num
                t_points(t_n) = t_n*dt;

                %Get the initial values of T and psi for this timestep
                T_n = permute(T(t_n, :, :), [2,3,1]);
                psi_nMinusHalf = permute(psi(t_n, :, :), [2,3,1]);

                %1. Calculate psi^n from Poisson equation (using SOR)      
                [dT_dz, dT_dx] = gradient(T_n, dx, dz);
                rhs_n = Rm * dT_dx;
                
                
                psi_n = poissonSolverTwoDimensions(rhs_n, psi_nMinusHalf, x_points, z_points, SOR_parameter);
                [u_n, v_n] = velocityFromPsi(psi_n, dx, dz);
                [psi_z, psi_x] = gradient(psi_n, dx, dz); 

                mod_u_2(rm_i, t_n) = nanmean(nanmean(abs(u_n).*abs(u_n) + abs(v_n).*abs(v_n)));

                %2. Calculate T^(n+1/2) from heat equation with a time step dt/2 using u^n  (ADI)
                %T_nPlusHalf = heatSolver(T_n, u_n, v_n, dt/2, x_points, z_points, sigma, k, true);
                T_nPlusHalf = heatSolverADI(x_points, z_points, x_num-1, z_num-1, dx, dz, dt/2, T_n, psi_z, psi_x, cartesian, neumann, T_x, T_z, T_xx, T_zz, T_xf_x);

                
                
                [dT_dz, dT_dx] = gradient(T_nPlusHalf, dx, dz);
                rhs_nPlusHalf = Rm * dT_dx;

                psi_nPlusHalf = poissonSolverTwoDimensions(rhs_nPlusHalf, psi_n, x_points, z_points, SOR_parameter);
                [u_nPlusHalf, v_nPlusHalf] = velocityFromPsi(psi_nPlusHalf, dx, dz);
                [psi_z, psi_x] = gradient(psi_n, dx, dz); 

                %4. Calculate T^(n+1) from heat equation with a time step dt using u^(n+1/2) for the velocity.
                if mod(t_n, 2) == 0
                    %fprintf('Calculating T^%d out of T^%d required iterations \n', t_n + 1, t_num);
                end
                %T_nPlusOne = heatSolver(T_n, u_nPlusHalf, v_nPlusHalf, dt, x_points, z_points, sigma, k, true);
                T_nPlusOne = heatSolverADI(x_points, z_points, x_num-1, z_num-1, dx, dz, dt/2, T_nPlusHalf, psi_z, psi_x, cartesian, neumann, T_x, T_z, T_xx, T_zz, T_xf_x);

                T(t_n + 1, :, :) = T_nPlusOne;
                psi(t_n + 1, :, :) = psi_nPlusHalf;
                
                %Calculate nusselt number
                [dT_dz, dT_dx] = gradient(T_nPlusOne, dx, dz);
                integrand = dT_dz;
                integral = trapz(x_points,integrand);
                Nus(t_n) = -integral(:, 1);

                %Check for steady state (compare T_nPlusOne with T_n)
                old_diff = diff;
                diff = max(max(abs(T_nPlusOne-T_n)));
                %if (diff > 1.5*old_diff)
                %    fprintf('Solution is diverging\n');
                %    break;
                %end

                frac_diff = abs(T_nPlusOne - T_n) ./ T_nPlusOne;
                max_frac_diff = max(max(frac_diff));                
                
                
                
                if (max_frac_diff < steady_state_criteria)
                    %If we've reached the steady state, leave the time loop
                    fprintf('Reached steady state\n');
                    T_final =  permute(T(t_n, :, :), [2,3,1]);      psi_final = permute(psi(t_n, :, :), [2,3,1]);
                    u_final = u_nPlusHalf;                          v_final = v_nPlusHalf;
                    break;
                else
                    
                    if (mod(t_n, 2) == 0)
                        %fprintf('Fractional difference with previous timestep: %1.10f \n', max_frac_diff);

                        if show_intermediate_graph
                            clf;
                            mesh(z_points, x_points, T_nPlusOne);
                            xlabel('z');
                            ylabel('x');
                            zlabel('T');
                            drawnow;
                        end
                    end

                    if t_n == t_num
                        T_final =  permute(T(t_n, :, :), [2,3,1]);      psi_final = permute(psi(t_n, :, :), [2,3,1]);
                        u_final = u_nPlusHalf;                          v_final = v_nPlusHalf;
                    end
                end

            end

        end
        
        T_frac_diff = (T_nPlusOne - T_analytic) ./ (T_analytic + 1);
        T_av_frac_diffs(x_num_i) = nanmean(nanmean(abs(T_frac_diff)));
        T_max_frac_diffs(x_num_i) = max(max(abs(T_frac_diff)));
        
        psi_frac_diff = abs(psi_nPlusHalf - psi_analytic) ./ (psi_analytic + 1);
        psi_av_frac_diffs(x_num_i) = nanmean(nanmean(psi_frac_diff));
        
        %mesh(T_analytic);
        %pause;
        %mesh(psi_frac_diff);
        %pause;
        %mesh(T_frac_diff);
        %pause;
        
    end
    
    mod_u_2_dt(dt_i) = nanmean(nanmean(abs(u_final).*abs(u_final) + abs(v_final).*abs(v_final)));
end


%Calculate the nusselt number
%Sum d(theta)/dz at z=0 over all x
[T_z, T_x] = gradient(T_final, dx, dz);
[psi_z, psi_x] = gradient(psi_final, dx, dz);

T = T_final;
v_z = psi_x;

%integrand = T_z + v_z .* T;
integrand = T_z;
integrand = (-3 * T(:, 1) + 4*T(:, 2) - T(:, 3))/(2*dz);
integral = trapz(x_points,integrand);
Nu = -integral(:, 1);

Nu_matlab = -trapz(x_points,T_z(:, 1));
%Nu_numeric = -sum(((-3 * T(:, 1) + 4*T(:, 2) - T(:, 3))/(2*dz)) * dx);
Nu_numeric = -trapz(x_points, ((-3 * T(:, 1) + 4*T(:, 2) - T(:, 3))/(2*dz)));
Nu_middle = -trapz(x_points,T_z(:, 15) + v_z(:,15).* T(:, 15));
Nu_test = -simpsum(T_z(:, 1), x_points);

%Nu = sum(T_z(:, 1))*dx/L;
fprintf('Ra: %3.2f \n', Rm);
fprintf('Nusselt number: %1.5f (matlab), %1.5f (numeric), %1.5f (middle)\n', Nu_matlab, Nu_numeric, Nu_middle);
fprintf('Max psi: %1.5f \n', max(max(psi_final)));

plot(t_points, Nus);
pause;


title_string = sprintf('T numeric (steady state). Rm: %3.1f, L: %2.2f, Nu: %1.3f', Rm, L, Nu);

%myPlot(x_points, z_points, T_final, 'T', title_string);

%plot the streamlines
%streamslice(z_points, x_points, v_final, u_final)

if numel(dts) == 1 && t_num > 1   
    figure;
    hold on;
    colour = ['--g' '--r' '--m' '--y' '--b'];

    tit_str  = '|u|^2 for Ra/Ra_c = ';
    for i = 1:numel(Rms)
        tit_str = strcat(tit_str, num2str(Rms(i)), ', ');
    end
    tit_str = strcat(tit_str, ' \epsilon = ', num2str(epsilon), ', \Delta t = ', num2str(dt));

    title(tit_str);
    xlabel('Time (s)');
    ylabel('|u|^2');
    for i = 1:numel(Rms)   
        % Get the time values for this dataset
        t = t_points(1:numel(mod_u_2(i, :)));
   
        plot(t, mod_u_2(i, :), colour(i));
    end

    %{
    for k = 1:numel(mod_u_2);
        fprintf('%1.9f ', mod_u_2(k));
        if mod(k, 10) == 0
            fprintf('...\n');
        end
    end
    %}

    pause;

end

if numel(dts) > 1
    plot(dts, mod_u_2_dt, '--*');
    xlabel('\Delta t');
    ylabel('|u|^2');
    tit_str = strcat('Value of |u(t=', num2str(t_max), ')|^2, Rm/Rm_c = ', num2str(Rms(end)));
    title(tit_str);
    
    pause;
end

%{
fprintf('dts = [');
for q = 1:numel(dts)
    fprintf(num2str(dts(q)));
end
fprintf('] \n mod_u_2_dt = [');
for q = 1:numel(mod_u_2_dt)
    fprintf(num2str(mod_u_2_dt(q)));
end
fprintf('] \n');
%}

%Compare with analytic solution
if t_num == 1
    plot(dxs.*dxs, T_av_frac_diffs, '--*');
    xlabel('(\Delta x)^2');
    ylabel('Fractional difference')
    title('Average fractional difference in T for one timestep');
    pause;
    plot(dxs.*dxs, psi_av_frac_diffs, '--*');
    xlabel('(\Delta x)^2');
    ylabel('Average frac diff')
    title('Average fractional difference in \psi for half a timestep');
    pause;
end

end

%{
myPlot(x_points, y_points, psi_analytic, '\psi', '\psi analytic (steady state)');
myPlot(x_points, y_points, psi_final, '\psi', '\psi numeric (steady state)');
myPlot(x_points, y_points, psi_diff_percent, '100 * \Delta \psi / \psi_{analytic}', strcat('\psi percentage difference', vars));

myPlot(x_points, y_points, T_analytic, 'T', 'T analytic (steady state)');

myPlot(x_points, y_points, T_diff_percent, '100 * \Delta T / T_{analytic}', strcat('T percentage difference', vars));
%}

%myPlot(x_points, y_points, u_analytic, 'u', 'u analytic (steady state)');
%myPlot(x_points, y_points, u_final, 'u', 'u numeric (steady state)');
%myPlot(x_points, y_points, u_diff, '\Delta u', 'u difference');

%myPlot(x_points, y_points, v_analytic, 'v', 'v analytic (steady state)');
%myPlot(x_points, y_points, v_final, 'v', 'v numeric (steady state)');
%myPlot(x_points, y_points, v_diff, '\Delta v', 'v difference');


 
%==========================================================================
% Plotting function to save a few lines
%==========================================================================

function myPlot(x_points, y_points, z,  z_title, title_text)
    %Figure1=figure(1);clf;
    %set(Figure1,'defaulttextinterpreter','latex');
    mesh(y_points, x_points, z);
    title(title_text);
    xlabel('z');
    ylabel('x');
    zlabel(z_title);
    pause;
end

%=========================================================================
% Calculate velocities from psi
%=========================================================================
function [u, v] = velocityFromPsi(psi_n, dx, dz)
    %Use inbuilt matlab function
    [psi_z, psi_x] = gradient(psi_n, dx, dz); %This IS correct (even though dPsiDx and dPsiDy seem to be the wrong way around)
    u = -psi_z;
    v = psi_x; 
    
    %We get spurious results on the edges, but these values aren't used in
    %calculations so not an issue.
end




