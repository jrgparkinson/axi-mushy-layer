%Plot over a continuous range of R, and a few discrete Rm values

clear; close all;
%H = 0.25;
theta_inf = 1.4; %Need to decide what this will be
R_vals = 0.001:0.003:0.4; %0.001:0.003:0.4
R_num = numel(R_vals);

r_num = 40;

lineColor = ['r', 'b', 'k'];

Rm_vals = [55 60 65]; %55:0.5:57;

Rm_num = numel(Rm_vals);

zero_line = nan * ones(R_num, Rm_num);

%flst=dir([pwd '/data/theta_inf_',theta_inf,'*.mat']);
%files_to_search={flst.name};

%Post processing
%Need two more lines to plot - one for when stuff is zero, the other to
%connect up to the bifurcation curve


for Rm_i = Rm_num:-1:1
    Rm = Rm_vals(Rm_i);
    
    flux(:, Rm_i) = nan * ones(R_num, 1);
    flux_unscaled(:, Rm_i) = nan * ones(R_num, 1);
    a_vals(:, Rm_i) = nan * ones(R_num, 1);
    
    for R_i = 1:numel(R_vals)
        R = R_vals(R_i);
        
        connecting_line(R_i, :) = [0 0];
        
        fprintf('----------------\n');

        
        %Need to find the file with this Rm and R and also the correct
        %theta_inf
        %See if we have already found the right file
        filelocation = file_name_with_vars3(R, Rm, theta_inf);
        
        if filelocation == false
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
            flux_unscaled(R_i, Rm_i) = nan;
            variation(R_i, Rm_i) = nan; 
        else
            load(filelocation);
            filename = strrep(filelocation, 'C:\Users\Jamie\Documents\MATLAB/data/theta_inf_1.4/','');
            fprintf('Loaded %s\n', filename);
            
            a_vals(R_i, Rm_i) = a(1);
            flux(R_i, Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1), C, 1);
            flux_unscaled(R_i, Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1), C, 0);
            deviation = abs(abs(theta - (-1 + z/constants('H')))./(theta-0.01));           
            variation(R_i, Rm_i) = mean(deviation(:, r_num));
            %fprintf('%1.5f, ',variation(R_i, Rm_i)); 
        end
        %{
        %If not, try and find it
        H = find_H(theta_inf, R, Rm, r_num, files_to_search);
        %H=nan; %Stop searching now to save time
        
        if isnan(H)
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
        else
            H_for_R(R) = H;
            
            filename = find_file_with_vars(R,H,Rm, r_num, files_to_search);
            
                
                load([pwd filename]);
            
                a_vals(R_i, Rm_i) = a(1);
                flux(R_i, Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));
          
        end
        %}
        
    end
    
    a_vals_Rm = a_vals(:, Rm_i);
    
    for R_i = (R_num-1):-1:1
        if ( a_vals_Rm(R_i) == 0 || isnan(a_vals_Rm(R_i))   ) && ...
                ( a_vals_Rm(R_i+1) > 0 )
            
            connecting_line(Rm_i, :) = [R_vals(R_i+1), a_vals_Rm(R_i+1)];
            
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
            flux_unscaled(R_i, Rm_i) = nan;
            
        end
    end

    for R_i = 1:(R_num-1)

        R = R_vals(R_i);
        if R <= connecting_line(Rm_i, 1)
            zero_line(R_i, Rm_i) = 0;
        else
            zero_line(R_i, Rm_i) = nan;
        end
    end
    
end

fig1 = figure(1); hold on;
x_vals = R_vals.^2;

for Rm_i = 1:numel(Rm_vals)
    
    a_plot = a_vals(:, Rm_i);
    flux_plot = flux(:, Rm_i);
    
    flux_unscaled_plot = flux_unscaled(:, Rm_i);
    
    variation_plot = variation(:, Rm_i);
    flux_plot = smooth(x_vals,flux_plot,0.04,'moving');
    
    idx = ~any(isnan(flux_plot.'),1);
    [hAx,hLine1(Rm_i),hLine2(Rm_i)] = plotyy(x_vals(idx), flux_unscaled_plot(idx), x_vals(idx), flux_plot(idx)); %flux_unscaled_plot
    
    xlabel('R^2')
    ylabel(hAx(1),'F') % left y-axis
    ylabel(hAx(2),'F/R') % right y-axis
    
    hLine1(Rm_i).LineStyle = '--';
    hLine2(Rm_i).LineStyle = '-';
    %hLine2(Rm_i).Marker = 'x';
    
    hLine1(Rm_i).Color=lineColor(Rm_i);
    hLine2(Rm_i).Color=lineColor(Rm_i);
    
    hLine2(Rm_i).DisplayName = strcat('Rm = ', num2str(Rm_vals(Rm_i)));
    hLine1(Rm_i).DisplayName = strcat('Rm = ', num2str(Rm_vals(Rm_i)));
    
    %Plot extra odd bits
    %{
    hLine3(Rm_i) = plot(R_vals, zero_line(:, Rm_i));
    hLine3(Rm_i).LineStyle = '-';
    hLine3(Rm_i).Color=lineColor(Rm_i);
    
    %Plot connecting line
    connecting_line_pos = connecting_line(Rm_i, 1);
    connecting_line_height = connecting_line(Rm_i, 2);
    
    connector = line([connecting_line_pos connecting_line_pos],[0 connecting_line_height],'Color',lineColor(Rm_i));
    connector.LineStyle = ':';
    %}
   
    axisLimits_a = [0 x_vals(end)*1.2 min(min(flux_unscaled_plot))*0.9 max(max(flux_unscaled_plot))*1.1];
    axisLimits_flux = [0 x_vals(end)*1.2 min(min(flux))*0.9 max(max(flux))*1.1];
    
    axis(hAx(1), axisLimits_a);
    axis(hAx(2), axisLimits_flux);
    
    %axis(hLine4(R_i), axisLimits_a);
    %axis(hLine5(R_i), axisLimits_flux);
    
    
end

l = legend([fliplr(hLine2)]);
l.Location = 'southeast';

%axis([40 70 0 0.06]);


%{
pause;

plot(Rm_vals, theta_inf);
xlabel('Rm'); ylabel('\theta_\infty');
axis([40 70 0 4]);
pause;
  

plot(Rm_vals, a_vals);
xlabel('Rm'); ylabel('a');
axis([40 70 0 0.04]);
%}

