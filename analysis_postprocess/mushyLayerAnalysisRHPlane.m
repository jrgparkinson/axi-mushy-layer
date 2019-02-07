%Results analysis
clear; close all;

Rm = 60;

[R_vals, H_vals] = meshgrid(0.19:0.006:0.45, 0.19:0.003:0.45);
%first index: Rm
%second index: R
%third index: H

R_num = numel(R_vals(1, :));
H_num = numel(H_vals(:, 1));


% Get all the data we have already processed
parsed_filename = strcat('\data\parsed_data_Rm',num2str(Rm),'.mat');
load([pwd parsed_filename]);
% This should return R_vals_p, H_vals_p, Rm_vals_p, psi_max_p,
psi_max_p = interp2(R_vals_p, H_vals_p,    psi_max_p,    R_vals, H_vals);

%a_vals_p = nan*psi_max_p;
%theta_inf_p = nan*psi_max_p;
%flux_p = nan*psi_max_p;

a_vals_p = interp2(R_vals_p, H_vals_p,    a_vals_p,    R_vals, H_vals);
theta_inf_p = interp2(R_vals_p, H_vals_p,    theta_inf_p,    R_vals, H_vals);
flux_p = interp2(R_vals_p, H_vals_p,    flux_p,    R_vals, H_vals);


flst=dir([pwd '/data/*.mat']);
files_to_search={flst.name};

%Get all the available data files so we can search them
%flst=dir([pwd '/data/*.mat']);
%flst={flst.name};

for R_i = 1:R_num
    R = R_vals(1, R_i);
    
    for H_i = 1:H_num
        H = H_vals(H_i, 1);
        
        
        %First check if we already have data
        if ~isnan(theta_inf_p(H_i, R_i))
            psi_max(H_i, R_i) = psi_max_p(H_i, R_i);            
            a_vals(H_i, R_i) = a_vals_p(H_i, R_i);
            theta_inf(H_i, R_i) = theta_inf_p(H_i, R_i);
            flux(H_i, R_i) = flux_p(H_i, R_i);
            
        else
            %if we don't have data, generate it
            
            filename = find_file_with_vars(R,H,Rm, files_to_search);
            
            if filename == false
                
                psi_max(H_i, R_i) = nan;
                a_vals(H_i, R_i) = nan;
                theta_inf(H_i, R_i) = nan;
                flux(H_i, R_i) = nan;
            else
                load([pwd filename]);
                fprintf('Loaded %s \n', filename);
                
                %Analyse the data
                psi_max(H_i, R_i) = max(max(psi));
                
                a_vals(H_i, R_i) = a(1);
                theta_inf(H_i, R_i) = calculateThetaInfinity(r,z,theta, psi);
                flux(H_i, R_i) = calculateSolidHeatFlux(r,z,theta,constants);
                
                %psi_max_Rm(Rm_i) = psi_max;
            end
        end
        
    end
end

% Save data
R_vals_p = R_vals;  H_vals_p = H_vals; 
psi_max_p = psi_max; theta_inf_p = theta_inf; a_vals_p = a_vals; flux_p = flux;
save([pwd parsed_filename],'R_vals_p','H_vals_p','psi_max_p','theta_inf_p','a_vals_p','flux_p');
fprintf('Saved data as %s \n', parsed_filename);


if R_num > 3 && H_num > 3
    [R_interp, H_interp] = meshgrid(0.22:0.005:0.37, 0.22:0.005:0.37);
    %psi_max_interp = interp2(R_vals(~isnan(psi_max)), H_vals(~isnan(psi_max)), ...
    %    psi_max(~isnan(psi_max)), R_interp, H_interp, 'spline');
    
    
    R_plot = R_vals;    H_plot = H_vals;    psi_max_plot = psi_max;
    
    mesh( R_plot, H_plot, psi_max_plot);
    xlabel('R'); ylabel('H'); 
    %title(strcat('\psi_{max} for Rm = ', num2str(Rm) ));
    pause;
    
    
    contour( R_plot, H_plot, psi_max_plot, 20);
    xlabel('R'); ylabel('H'); 
    %title(strcat('\psi_{max} for Rm = ', num2str(Rm) ));
    pause;
    
    contour( R_plot, H_plot, a_vals, 20);
    xlabel('R'); ylabel('H'); 
    pause;
    
    contour( R_plot, H_plot, theta_inf, 20);
    xlabel('R'); ylabel('H'); 
    pause;
    
    contour( R_plot, H_plot, flux, 20);
    xlabel('R'); ylabel('H'); 
    pause;
    
    
end





