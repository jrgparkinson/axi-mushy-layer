%Plot over a continuous range of Rm, and a few discrete R values

% R = 0.22: lowest Rm = 54
% R = 0.28: lowest Rm = need to investigate
% R = 0.3: lowest Rm =

clear; close all;
H = 0.25;
theta_inf = 1.4; %Need to decide what this will be
Rs = [0.22 0.26 0.3];  %0.28 0.34
r_num = 40;

lineColor = ['r', 'b', 'k', 'm', 'y'];

Rm_vals = [45:0.5:75]; %55:0.5:57;

Rm_num = numel(Rm_vals);

flst=dir([pwd '/data/*.mat']);
files_to_search={flst.name};

%Post processing
%Need two more lines to plot - one for when stuff is zero, the other to
%connect up to the bifurcation curve

for R_i = 1:numel(Rs)
    R = Rs(R_i);
    
    connecting_line(R_i, :) = [0 0];
    
    %load(strcat('H_for_Rm_map_R',num2str(R),'.mat'));   
    H_for_Rm = containers.Map({0}, [0]);

    %H_for_Rm(56) = 0.29210;
    
    for Rm_i = Rm_num:-1:1
        Rm = Rm_vals(Rm_i);
        
        fprintf('----------------\n');
        
        flux(R_i, Rm_i) = nan;
        a_vals(R_i, Rm_i) = nan;
        
        %Need to find the file with this Rm and R and also the correct
        %theta_inf
        filelocation = file_name_with_vars3(R, Rm, theta_inf);
        
        if filelocation == false
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
        else
            load(filelocation);
            filename = strrep(filelocation, 'C:\Users\Jamie\Documents\MATLAB/data/theta_inf_1.4/','');
            fprintf('Loaded %s\n', filename);
            
            a_vals(R_i, Rm_i) = a(1);
            flux(R_i, Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1),C, 1);
        end
        
        %{
        if isnan(H)
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
        else
            H_for_Rm(Rm) = H;
            
            filename = find_file_with_vars(R,H,Rm, r_num, files_to_search);
            load([pwd filename]);
            
            a_vals(R_i, Rm_i) = a(1);
            flux(R_i, Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));
        end
        %}
         
        if Rm_i < numel(a_vals(R_i, :)) && ...
                ( a_vals(R_i, Rm_i) == 0 || isnan(a_vals(R_i, Rm_i))   ) && ...
                ( a_vals(R_i, Rm_i+1) > 0 )
            
            connecting_line(R_i, :) = [Rm_vals(Rm_i+1), a_vals(R_i, Rm_i+1)];
            
            a_vals(R_i, Rm_i) = nan;
            flux(R_i, Rm_i) = nan;
            
        end
    end
    
    for Rm_i = 1:Rm_num
        Rm = Rm_vals(Rm_i);
        if Rm <= connecting_line(R_i, 1) 
            zero_line(R_i, Rm_i) = 0;
        else
            zero_line(R_i, Rm_i) = nan;
        end
    end
    
    %Save H data
    save(strcat('H_for_Rm_map_R',num2str(R),'.mat'),'H_for_Rm');
end

fig1 = figure(1); hold on;

x_vals = Rm_vals;

for R_i = 1:numel(Rs)
    
    a_plot = a_vals(R_i, :);
    flux_plot = flux(R_i, :);
    
    

    %flux_plot = smooth(Rm_vals,flux_plot,0.05,'moving');
    
    idx = ~any(isnan(a_plot),1);
    [hAx,hLine1(R_i),hLine2(R_i)] = plotyy(x_vals(idx), a_plot(idx), x_vals(idx), flux_plot(idx));
    
    xlabel('Rm')
    ylabel(hAx(1),'a') % left y-axis
    ylabel(hAx(2),'F/R') % right y-axis
    
    hLine1(R_i).LineStyle = '--';
    %hLine1(R_i).Marker = '+';
    hLine2(R_i).LineStyle = '-';
    %hLine2(R_i).Marker = 'x';
    
    hLine1(R_i).Color=lineColor(R_i);
    hLine2(R_i).Color=lineColor(R_i);
    
    hLine2(R_i).DisplayName = strcat('R = ', num2str(Rs(R_i)));
    hLine1(R_i).DisplayName = strcat('R = ', num2str(Rs(R_i)));
    
    %Plot extra odd bits
    hLine3(R_i) = plot(Rm_vals, zero_line(R_i, :));
    hLine3(R_i).LineStyle = '-';
    hLine3(R_i).Color=lineColor(R_i);
   
    
    %Plot connecting line
    connecting_line_pos = connecting_line(R_i, 1);
    connecting_line_height = connecting_line(R_i, 2);
    
    connector = line([connecting_line_pos connecting_line_pos],[0 connecting_line_height],'Color',lineColor(R_i));
    connector.LineStyle = ':';
    
    axisLimits_a = [x_vals(1)*0.9 x_vals(end)*1.2 0 0.04];
    axisLimits_flux = [x_vals(1)*0.9 x_vals(end)*1.2 0 max(max(flux))*1.2];
    
    axis(hAx(1), axisLimits_a);
    axis(hAx(2), axisLimits_flux);
    
    %axis(hLine4(R_i), axisLimits_a);
    %axis(hLine5(R_i), axisLimits_flux);
    
    
end

l = legend([hLine2]);
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

