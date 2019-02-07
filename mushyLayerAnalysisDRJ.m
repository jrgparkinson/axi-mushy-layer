%Plot over a continuous range of R, and a few discrete Rm values

clear; close all;
theta_inf = 1.4;
R_vals = [0.22 0.26 0.3];  R_num = numel(R_vals);
r_num = 40;
Rm_vals = 45:0.5:75;
Rm_num = numel(Rm_vals);

r_R = r_num;

for R_i = 1:R_num
    R = R_vals(R_i);
    
    for Rm_i = 1:Rm_num
        Rm = Rm_vals(Rm_i);
        
        filelocation = file_name_with_vars3(R, Rm, theta_inf);
        
        
        if filelocation == false

            variation(R_i, Rm_i, :) = nan*[1:r_num];
            %fprintf('Could not find %s \n', filelocation);

        else
            load(filelocation);
            %fprintf('Loaded %s \n', filename);
            
            %Analyse the data
            variation(R_i, Rm_i, :) = calculateNonlinearTheta(r, z, constants('H'), theta, 1:r_num); 
        end
        
        %Rm_R_mean = mean(variation(R_i, Rm_i, :)); 
        %fprintf('Rm = %2.1f, R=%1.5f, average deviation = %1.5f\n', Rm, R, Rm_R_mean);

    end
    
    
end

fig1 = figure(1); hold on;

r_plots = [40];

for R_i = 1:R_num
    
    for r_plot_i = 1:numel(r_plots)
        r_plot = r_plots(r_plot_i);
        
        %variation_plot = variation(R_i, :,r_plot);
        variation_plot = variation(R_i, :, r_plot);
        
        hLines(R_i, r_plot_i) = plot(Rm_vals, variation_plot,'x');
        xlabel('Rm');
        ylabel('\delta');
        hLines(R_i, r_plot_i).DisplayName = strcat('R = ', num2str(R_vals(R_i)));
        
        axis([40 80 0 max(variation_plot)*1.3]);
        box on;
        
    end
end

l = legend([hLines]);
l.Location = 'southeast';
