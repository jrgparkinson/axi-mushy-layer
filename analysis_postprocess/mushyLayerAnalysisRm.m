clear; close all;
H = 0.25;
R = 0.28; 
R2 = 0.31;
R3 = 0.34;
r_num = 40;

Rm_vals = 40:0.5:80;

Rm_num = numel(Rm_vals);

flst=dir([pwd '/data/*.mat']);
files_to_search={flst.name};


for Rm_i = 1:Rm_num
    Rm = Rm_vals(Rm_i);
 
    filename = find_file_with_vars(R,H,Rm, r_num, files_to_search);
    
    if filename == false 
        %psi_max(Rm_i) = nan; 
        %theta_inf(Rm_i) = nan;
        flux1(Rm_i) = nan;
        a_vals1(Rm_i) = nan;
    else
        load([pwd filename]);
        fprintf('Loaded %s \n', filename);
        
        %Analyse the data
        %psi_max(Rm_i) = max(max(psi));
        %theta_inf(Rm_i) = calculateThetaInfinity(r,z,theta, psi);
        a_vals1(Rm_i) = a(1);
        flux1(Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));

        %psi_max_Rm(Rm_i) = psi_max;
    end
    
    %Get data for another R value
    filename2 = find_file_with_vars(R2,H,Rm, r_num, files_to_search);
    
    if filename2 == false 
        %psi_max(Rm_i) = nan; 
        %theta_inf(Rm_i) = nan;
        flux2(Rm_i) = nan;
        a_vals2(Rm_i) = nan;
    else
        load([pwd filename2]);
        fprintf('Loaded %s \n', filename2);
        
        %Analyse the data
        %psi_max(Rm_i) = max(max(psi));
        %theta_inf(Rm_i) = calculateThetaInfinity(r,z,theta, psi);
        a_vals2(Rm_i) = a(1);
        flux2(Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));

        %psi_max_Rm(Rm_i) = psi_max;
    end
    
    %Get data for another R value
    filename3 = find_file_with_vars(R3,H,Rm, r_num, files_to_search);    
    if filename3 == false 
        flux3(Rm_i) = nan;
        a_vals3(Rm_i) = nan;
    else
        load([pwd filename3]);
        fprintf('Loaded %s \n', filename3);
        a_vals3(Rm_i) = a(1);
        flux3(Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));
    end
    
    
end

%Post processing
for Rm_i = (round(Rm_num/2)):-1:1
    if a_vals1(Rm_i) == 0 || isnan(a_vals1(Rm_i))
        a_vals1(Rm_i) = 0;        
        flux1(Rm_i) = 0;
        
    end
    
    if a_vals2(Rm_i) == 0  || isnan(a_vals2(Rm_i))
        a_vals2(Rm_i) = 0;        
        flux2(Rm_i) = 0;
        
    end
    
    if a_vals3(Rm_i) == 0 || isnan(a_vals3(Rm_i))
        a_vals3(Rm_i) = 0;        
        flux3(Rm_i) = 0;
    end
end
    
fig1 = figure(1); hold on;
%axis([40 70 0 1]);

%[hAx,hLine1,hLine2] = plotyy(Rm_vals, a_vals, Rm_vals, flux);
%[hAx,hLine1,hLine2] = plotyy([Rm_vals, Rm_vals], [a_vals1, a_vals2], ...
%    [Rm_vals', Rm_vals'], [flux1', flux2']);
[hAx,hLine1,hLine2] = plotyy(Rm_vals, a_vals1, Rm_vals, flux1);
[hAx_2,hLine1_2,hLine2_2] = plotyy(Rm_vals, a_vals2, Rm_vals, flux2);
[hAx_3,hLine1_3,hLine2_3] = plotyy(Rm_vals, a_vals3, Rm_vals, flux3);

xlabel('Rm')
ylabel(hAx(1),'a') % left y-axis
ylabel(hAx(2),'Flux') % right y-axis

hLine1.LineStyle='--';
hLine1_2.LineStyle='--';
hLine1_3.LineStyle='--';

hLine2.LineStyle = '-';
hLine2_2.LineStyle = '-';
hLine2_3.LineStyle = '-';

hLine1.Color='red';
hLine2.Color='red';

hLine1_2.Color='blue';
hLine2_2.Color='blue';

hLine1_3.Color='black';
hLine2_3.Color='black';

axisLimits_a = [40 75 0 0.04];
axisLimits_flux = [40 75 0 30];
axis(hAx(1), axisLimits_a);
axis(hAx(2), axisLimits_flux);

axis(hAx_2(1), axisLimits_a);
axis(hAx_2(2), axisLimits_flux);

axis(hAx_3(1), axisLimits_a);
axis(hAx_3(2), axisLimits_flux);

hLine1.DisplayName = 'R = 0.25';
hLine1_2.DisplayName = 'R = 0.31';
hLine1_3.DisplayName = 'R = 0.34';

legend([hLine1 hLine1_2 hLine1_3]);

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

