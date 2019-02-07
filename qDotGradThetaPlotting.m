%Plot q.grad(theta) for an example steady state and illustrate the 
%issues with linear extrapolation

close all;
clear;

load([pwd strcat('\data\mushyLayerPrevSteadyState_relaxed_Points40Rm60H0.271R0.274b0.034282.mat')]);
R = constants('R');
fprintf('R = %1.5f\n',R);

[dr, dz, r_num, z_num] = meshGridProperties(r,z);

[theta_r, theta_z] = gradient2order(theta, dr, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r, z, psi);
q_r = u_r;  q_z = u_z-1;

qDotGradTheta = theta_r .* q_r + theta_z.*q_z + 3.8;

z_i = round(0.5*z_num/3);
q_i = qDotGradTheta(z_i, :);
r_i = r(z_i, :);

b = a(1) + dr; %approx

dr_fine = 0.001;
r_fine = a:dr_fine:R;
q_fine = interp1(r_i, q_i, r_fine, 'spline');

figure;
hold on;
box on;
r_min = 0; r_max = R/2; r_range = r_max - r_min;
q_min = -0.25; q_max = 1.25; q_range = q_max-q_min;

largeAxis = [r_min r_max q_min q_max];
axis(largeAxis);
plot(r_fine, q_fine);
xlabel('r'); ylabel('q . \nabla \theta');

%Calculate some extrapolations
for j = 0:1
    i = 1+j*11;
   gradient = (q_fine(i+1) - q_fine(i))/dr_fine;
   r_points = 0:dr_fine:r_fine(i+1);
   extrap_line = (r_points-r_fine(i))*gradient + q_fine(i);
   plot(r_points, extrap_line, '--r');
   
   
   %Highlight ranges;
   x1 = r_fine(i);
   line([x1 x1],[-0.5 q_fine(i)],'Color',[0.7 0.7 0.7],'LineStyle','-');
   fprintf('extrapolation from r = %1.5f\n', x1);
   
   x0 = r_fine(i) - q_fine(i)/gradient;
   line([x0 x0],[-0.5 0],'Color',[0.5 0.5 0.5],'LineStyle','-');
   
   %Add an x axis
   r_points_full = 0:0.01:R;
   xaxis = r_points_full*0;
   plot(r_points_full, xaxis, 'k');
   
   
   
end

%insert full figure into top right corner

insert_r = 0.6;
insert_q = 0.6;
insert_width = 0.895 - insert_r;
insert_height = 0.915 - insert_q;

handaxes2 = axes('Position', [insert_r insert_q insert_width insert_height]);
plot(r_fine, q_fine, 'b')
rectangle('Position',[r_min q_min r_range q_range], 'LineStyle','--')
axis([0 R*1.01 -2.5 1.5]);
%set(handaxes2, 'Box', 'off')
%xlabel('r'); ylabel('q . \nabla \theta');




