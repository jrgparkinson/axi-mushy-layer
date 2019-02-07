%Plot over a continuous range of Rm, for the R value with max flux

clear; close all;
theta_inf = 1.4; %Need to decide what this will be
r_num = 40;

Rm_vals = [45:0.5:70]; %55:0.5:57;

Rm_num = numel(Rm_vals);

for Rm_i = Rm_num:-1:1
    Rm = Rm_vals(Rm_i);
    
    fprintf('Rm = %2.1f \n', Rm);
    
    flux(Rm_i) = nan;
    a_vals(Rm_i) = nan;
    
    %Need to find the file with this Rm, max flux and also the correct
    %theta_inf
    folder = [pwd '/data/theta_inf_' num2str(theta_inf) '/'];

    file = strcat('mushyLayerPrevSteadyState_relaxed_Points',num2str(r_num),'Rm',num2str(Rm),'H*R*b*.mat');
    filesearch = [folder file];
    filelist = dir(filesearch);
    
    
    if numel(filelist) > 0
        
        if numel(filelist) > 1
            %Find the file with the largest solute flux
            largest_flux = 0;
            i_largest = 1;
            
            for i = 1:numel(filelist)
                load([folder filelist(1).name]);
                %%%
                flux_file(i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1), C, 1);
                
                if flux_file(i) > largest_flux
                    largest_flux = flux_file(i);
                    i_largest = i;
                end
            end
   
            filelocation = [folder filelist(i_largest).name];
            
        else
            %Only one file found
            filelocation = [folder filelist(1).name];
        end
        
    else
        filelocation = false;
    end
 
    
    if filelocation == false
        a_vals(Rm_i) = nan;
        flux(Rm_i) = nan;
    else
        load(filelocation);
        filename = strrep(filelocation, 'C:\Users\Jamie\Documents\MATLAB/data/theta_inf_1.4/','');
        fprintf('Loaded %s\n', filename);
        
        a_vals(Rm_i) = a(1);
        flux( Rm_i) = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1),C, 1);
        theta_inf_file = calculateThetaInfinity(r,z,theta,psi(:,:));
        
        fprintf('theta_inf = %1.5f\n', theta_inf_file);
    end
  
end


%Save H data
%save(strcat('H_for_Rm_map_R',num2str(R),'.mat'),'H_for_Rm');

fig1 = figure(1); hold on;

x_vals = Rm_vals;

flux = smooth(x_vals,flux,0.1,'moving');

idx = ~any(isnan(a_vals),1);
[hAx,hLine1,hLine2] = plotyy(x_vals(idx), a_vals(idx), x_vals(idx), flux(idx));

xlabel('Rm')
ylabel(hAx(1),'a') % left y-axis
ylabel(hAx(2),'F/R') % right y-axis

hLine1.LineStyle = '--';
%hLine1.Marker = '+';
hLine2.LineStyle = '-';
hLine2.Marker = 'x';

axisLimits_a = [x_vals(1)*0.9 x_vals(end)*1.2 0 0.04];
axisLimits_flux = [x_vals(1)*0.9 x_vals(end)*1.2 0 max(max(flux))*1.2];

axis(hAx(1), axisLimits_a);
axis(hAx(2), axisLimits_flux);


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

