%Produce plots of a single steady state
function mushyLayerPlotSteadyState(Rm, R, H, points)
%clear;
close all;


flst=dir([pwd '/data/*.mat']);
flst={flst.name};
filename = find_file_with_vars(R, H, Rm, points, flst);

if filename ~= false
    load([pwd filename]);
    
    R = constants('R'); H = constants('H');
    
    %Extrapolate solution
    [r_ext, z_ext, theta_ext, psi_ext, Q, i_min, b_i] = extrapolateSolution(r, z, theta, psi(:,:), a, C, constants);
    psi_ext_q = psi_ext - 0.5*r_ext.^2;
    
    [dr, dz, r_num, z_num] = meshGridProperties(r_ext, z_ext);
    

    %mesh(r_ext,z_ext, psi_ext); pause;
    
    % Find where the chimney should actually be
    Q = calculateQ(r_ext, z_ext,theta_ext,psi_ext);
    a_pred = zeros(z_num);
    
    for z_i = 1:z_num
        a_pred(z_i) = nan;
        for r_i = 2:r_num
            if (Q(z_i, r_i) > 0 && Q(z_i, r_i-1) < 0) ||...
                 (Q(z_i, r_i) < 0 && Q(z_i, r_i-1) > 0)   
                %Found the root
                a_pred(z_i) = r_ext(1, r_i-1 ) + ...
                    (Q(z_i, r_i) - Q(z_i, r_i-1)) * (r_ext(1, r_i) - r_ext(1, r_i-1 )) ;
                break;
            end
        end 
    end
    %plot(a_pred, z_ext(:, 1));
    %xlabel('a'); ylabel('z');
    %pause;
    
 
    %mesh(r_ext, z_ext, psi_ext_q); pause;
    
    %Switch geometry
    theta = flipud(theta); 
    r_ext = flipud(r_ext);
    z_ext = flipud(z_ext);
    
    %Calculate velocities
    
    %What is the maximum velocity within the mushy layer?
    a_i = round(numel(r_ext(1, :))*a(1)/R)+2;
    [u_r, u_z] = axiVelocitiesFromPsi(r_ext(:, a_i:end), z_ext(:, a_i:end), psi_ext(:, a_i:end));
    fprintf('Maximum u_r: %1.5f, u_z: %1.5f\n',max(max(u_r)), max(max(u_z)));
    
    %What is the solute flux?
    Flux = calculateSolidHeatFlux(r,z,theta,psi,constants,a(1));
    fprintf('Solute flux = %1.5f\n',Flux);
    
    
    %[q_r, q_z] = axiVelocitiesFromPsi(r_ext, z_ext, psi_ext_q);
    
    startr = 0:0.002:R; startz = -H*ones(size(startr));

    z_ext = -flipud(z_ext);
   

    fig = figure(1); hold on; box on;
    %set(fig,'defaulttextinterpreter','latex');
    
    % theta contours
    %[C2, h2] = contourf(r_ext, z_ext, theta_ext, 100);
    %colormap(jet);
    %h2.LineStyle = 'none';
    
    [C, h] = contour(r_ext, z_ext, theta_ext, 9); 
    h.LineStyle = '--';
    h.LineColor = 'black';
    h.LineWidth = 1.2;
    
    %cbar = colorbar;
    %ylabel(cbar, '\theta');
    
    % psi streamlines
    %streamline(r_ext, z_ext, u_r, u_z, startr, startz);
    psi_to_plot = psi_ext;

    numContours = round((max(max(psi_to_plot)) - 0)/0.004);
    fprintf('Using %d contours\n', numContours);
    [C3, h3] = contour(r_ext, z_ext, psi_to_plot, numContours);
    h3.LineStyle = '-';
    h3.LineColor = 'blue';
    h3.LineWidth = 1.2;
    %h3.LevelStep = 0.001;
    %streamline(r_ext, z_ext, q_r, q_z, startr, startz);
    
    %chimney
    a_point=a(1);
    chimneyWall = line([a_point a_point],[0 -H],'Color','red');
    chimneyWall.LineWidth = 1.2;
    
    % axes
    xlabel('$r$','Interpreter','LaTex','FontSize',20); ylabel('$z$','Interpreter','LaTex','FontSize',20);
    
   

    
else
    fprintf('Could not find a file with those variables\n');
end