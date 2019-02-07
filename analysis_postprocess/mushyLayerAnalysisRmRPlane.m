%Results analysis
clear; close all;

H = 0.2;

[R_vals, Rm_vals] = meshgrid(0.16:0.006:0.352, 40:0.5:80);
%first index: Rm
%second index: R
%third index: H

R_num = numel(R_vals(1, :));
Rm_num = numel(Rm_vals(:, 1));


% Get all the data we have already processed
parsed_filename = strcat('\data\parsed_data_H',num2str(H),'.mat');
%load([pwd parsed_filename]);
%psi_max_p = interp2(R_vals_p, Rm_vals_p,    psi_max_p,    R_vals, Rm_vals);

psi_max_p = nan*R_vals;

flst=dir([pwd '/data/*.mat']);
files_to_search={flst.name};

%Get all the available data files so we can search them
%flst=dir([pwd '/data/*.mat']);
%flst={flst.name};

for R_i = 1:R_num
    R = R_vals(1, R_i);
    
    for Rm_i = 1:Rm_num
        Rm = Rm_vals(Rm_i, 1);
        
        
        %First check if we already have data
        if ~isnan(psi_max_p(Rm_i, R_i))
            psi_max(Rm_i, R_i) = psi_max_p(Rm_i, R_i);
        else
            %if we don't have data, generate it
            
            filename = find_file_with_vars(R,H,Rm, 40, files_to_search);
            
            if filename == false
                
                psi_max(Rm_i, R_i) = nan;
            else
                load([pwd filename]);
                fprintf('Loaded %s \n', filename);
                
                %Analyse the data
                psi_max(Rm_i, R_i) = max(max(psi));
                
                %psi_max_Rm(Rm_i) = psi_max;
            end
        end
        
    end
end

% Save data
R_vals_p = R_vals;  Rm_vals_p = Rm_vals; 
psi_max_p = psi_max;
save([pwd parsed_filename],'R_vals_p','Rm_vals_p','psi_max_p');
fprintf('Saved data as %s \n', parsed_filename);


if R_num > 3 && Rm_num > 3
    %[R_interp, Rm_interp] = meshgrid(0.22:0.005:0.37, 0.22:0.005:0.37);
    %psi_max_interp = interp2(R_vals, H_vals, psi_max(), R_interp, H_interp, 'linear');
    
    R_plot = R_vals;    Rm_plot = Rm_vals;   psi_max_plot = psi_max;
    
    mesh( R_plot, Rm_plot, psi_max_plot);
    xlabel('R'); ylabel('Rm'); 
    %title(strcat('\psi_{max} for H = ', num2str(H) ));
     c = colorbar();
    c.Label.String = '\psi_{max}'; 
    pause;
    
    
    contour( R_plot, Rm_plot, psi_max_plot, 20);
    xlabel('R'); ylabel('Rm'); 
    %title(strcat('\psi_{max} for H = ', num2str(H) ));
    c = colorbar();
    c.Label.String = '\psi_{max}';
    pause;
    
    
end




