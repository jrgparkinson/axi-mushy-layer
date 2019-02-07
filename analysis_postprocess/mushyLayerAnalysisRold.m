%Plot over a continuous range of R, and a few discrete Rm values

clear; close all;
H = 0.25; %0.2
r_num = 40;
Rm_vals = [60]; %[57 62 67 72]
Rm_num = numel(Rm_vals);
Rm_color = ['r', 'b', 'k', 'm'];

R_vals = 0.01:0.003:0.455;
R_num = numel(R_vals);

flst=dir([pwd '/data/*.mat']);
files_to_search={flst.name};

for Rm_i = 1:Rm_num
        Rm = Rm_vals(Rm_i);
        
        for R_i = 1:R_num
            R = R_vals(R_i);
   
        filename = find_file_with_vars(R,H,Rm, r_num, files_to_search);
        
        if filename == false
            
            psi_max(Rm_i, R_i) = nan;
            a_vals(Rm_i, R_i) = nan;
            theta_inf(Rm_i, R_i) = nan;
            flux(Rm_i, R_i) = nan;
            %fprintf('Could not find %s \n', filename);
            
        else
            load([pwd filename]);
            %fprintf('Loaded %s \n', filename);
            
            %Analyse the data
            %psi_max(Rm_i, R_i) = max(max(psi));
            a_vals(Rm_i, R_i) = a(1);
            theta_inf = calculateThetaInfinity(r,z,theta, psi);
            fprintf('Rm = %2.2f, R = %1.5f, H = %1.5f, theta_inf = %1.4f \n',Rm,R, H, theta_inf);
            flux(Rm_i, R_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));
            
            %psi_max_Rm(Rm_i) = psi_max;
            
            
            
        end
    end
    
end

%Post processing
for Rm_i = 1:Rm_num
    empty  = nan*[1:R_num];
    zero_line(Rm_i, :) = empty;
    
    
    for R_i = (R_num-1):-1:1
       
        if a_vals(Rm_i, R_i) == 0 || isnan(a_vals(Rm_i, R_i))
            a_vals(Rm_i, R_i) = nan;
            flux(Rm_i, R_i) = nan;
                
            
                
                
            if a_vals(Rm_i, R_i + 1) > 0
                zero_line(Rm_i, R_i) = 0;
                zero_line(Rm_i, R_i+1) = 0;
                connecting_line(Rm_i, :) = [R_vals(R_i+1), a_vals(Rm_i, R_i+1)]; 
                
                
            end
            
            if zero_line(Rm_i, R_i+1) == 0
                zero_line(Rm_i, R_i) = 0;
            end
                
            
            
        end
        
    end
end
    

fig1 = figure(1); hold on;

x_data = R_vals; % .^2

for Rm_i = 1:Rm_num

[hAx,hLine1(Rm_i),hLine2(Rm_i)] = plotyy(x_data, a_vals(Rm_i,:), x_data, flux(Rm_i,:));

xlabel('R^2')
ylabel(hAx(1),'a') % left y-axis
ylabel(hAx(2),'Flux') % right y-axis

hLine1(Rm_i).Color = Rm_color(Rm_i);
hLine2(Rm_i).Color = Rm_color(Rm_i);

hLine1(Rm_i).LineStyle = '--';

%Plot connecting line
connecting_line_pos = connecting_line(R_i, 1).^2;
connecting_line_height = connecting_line(R_i, 2);
connector = line([connecting_line_pos connecting_line_pos],[0 connecting_line_height],'Color',Rm_color(R_i));
connector.LineStyle = ':';

%Plot zero line
zLine = plot(x_data, zero_line(Rm_i, :));
zLine.Color = Rm_color(Rm_i);

hLine1(Rm_i).DisplayName = strcat('Rm = ', num2str(Rm_vals(Rm_i)));

axis(hAx(1), [0 1.2*x_data(end) 0 0.036]);
axis(hAx(2), [0 1.2*x_data(end) 0 40]);

set(hAx(1),'box','off');

%axis([40 70 0 0.06]);

end

l = legend([hLine1]);
l.Location = 'southeast';

