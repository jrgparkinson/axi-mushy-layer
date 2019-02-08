function [a_new, a_extrap] = extrapolate_a(r, z, t, a_n, constants, psi, theta, C_n)
[dr, dz, ~, z_num] = meshGridProperties(r, z);
[dt, ~] = gridProperties(t);

b = constants('b');
a_fixed = constants('a_fixed');
chimney_method = constants('chimney_method');

a_extrap = nan*ones(z_num, 1); %In case we don't use this var

%Set a_fixed < 0 to use relaxation
if a_fixed > 0
    a_new = a_fixed * ones(z_num, 1);
else
    %f extrap method
    if chimney_method == 1
        psi_q = psi - 0.5*r.^2;
        
        [theta_r, theta_z] = gradient2order(theta, dr, dz);
        %[u, v] = axiVelocitiesFromPsi(r, z, psi);
        %q_r = u;
        %q_z = v;
        [q_r, q_z] = axiVelocitiesFromPsi(r, z, psi_q);
        
        f_r = q_r.*theta_r;
        f_z = q_z.*theta_z;
        f_q = f_r + f_z;
        
        
        [u, v] = axiVelocitiesFromPsi(r, z, psi);
        f_u = u.*theta_r + v.*theta_z;
        
        %{
contour(r, z, theta.');
title('\theta');
ylabel('z');
xlabel('r');
pause;

mesh(r, z, theta_r.');
title('\theta');
ylabel('z');
xlabel('r');
pause;



mesh(r, z, f_q);
title('q. \nabla \theta');
xlabel('r');
ylabel('z');
pause;

mesh(r, z, f_u);
title('u. \nabla \theta');
xlabel('r');
ylabel('z');
pause;

f_mid = f_q( round(z_num/2) , :);
f_bottom = f_q( round(z_num) , :);
f_end = f_q( round(1) , :);
plot(r(1, :), f_mid, 'b', r(1, :), f_end, 'g', r(1, :), f_bottom, 'r');
xlabel('r')
ylabel('q. \nabla \theta');
title('q. \nabla \theta at z = 0 (red), H/2 (blue), H (green)');
pause;
        %}
        
        %df = (-3*f(:, 1) + 4*f(:, 2) - f(:, 3)) / (2*dr);
        [df, ~] = gradient2order(f_q, dr, dz);
        
        a_extrap = b - f_q(:, 1)./df(:, 1);
        
        %Q contour method
    elseif chimney_method == 2
        
        [r_concat, z_concat, T_concat, psi_concat, Q, i_min, b_i] = extrapolateSolution(r, z, theta, psi, a_n, C_n, constants);
        
        %mesh(r_concat(1:end-5, :), z_concat(1:end-5, :), Q(1:end-5, :)); pause;
        %plot(r_concat(10, :), Q(10, :)); pause;
        
        %Manually find where Q=0;
        a_interp = zeros(z_num, 1);
        for j = 1:z_num
            
            Q_j = Q(j, i_min+1:b_i);
            R = r_concat(1, i_min+1:b_i);
            a_interp(j, 1) = interp1(Q_j, R, 0);
            %fprintf('a(j = %d) = %1.2f \n', j, a_interp);
            %plot(r_concat(j, i_min+1:i_max+5), Q(j, i_min+1:i_max+5)); pause;
            
            
            
            clf;
            figure
            mesh(r_concat(1:end-5, :), z_concat(1:end-5, :), Q(1:end-5, :))
            xlabel('r'); ylabel('z'); title(strcat('Q = q.\nabla\theta steady state, last 5 z-points omitted, ', vars));
            savefig(strcat('Q',filename_vars, '.fig'));
            %pause;
            
            
        end
        
        %plot(a_interp, z(:, 1)); pause;
        
        
    %cfit method
    elseif chimney_method == 3 

        %Calculate q.grad(theta) at z=2H/3
        j = round(0.66*z_num);

        f_q = calculateQ(r,z,theta,psi);
   
        q = f_q(j, 1:4).';
        r_4 = r(j, 1:4).';
        
        % Fit to a cubic curve
        fo = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[10 0.006 0.05]);
           
        f = fittype('a*x^2+b*x+c','options', fo);
        fitobject = fit(r_4, q, f);
        coeffs = coeffvalues(fitobject);
            
        %a1 = coeffs(1); a2=coeffs(2); a3=coeffs(3); a4 = coeffs(4);
        a1=coeffs(1); a2=coeffs(2); a3 = coeffs(3);
        
        % Plot curve over 0 < r < 1.5*b on a finely sampled grid
        dr_extrap = 0.0001;
        r_extrap = 0:dr_extrap:(b*1.5); 
        q_extrap = (a1*r_extrap.^2 + a2*r_extrap + a3);

        
        %Find roots by change of sign
        b_i = round(b/dr_extrap) + 1;
        
        root = 0;
        for i = b_i:-1:2
            if q_extrap(i) > 0 && q_extrap(i-1) < 0
                % Interpolate between the two r points either side of the
                % change of sign
                root = r_extrap(i) - dr_extrap*q_extrap(i)/(q_extrap(i) - q_extrap(i-1));
                break;
            end            
        end
        
        a_extrap = root* ones(z_num, 1);
        %fprintf('extrapolated a: %1.6f \n', root);
        
        %{
figure(1); hold on;
        plot(r(j, :), f_q(j, :), 'b');
        plot(fitobject, r_4, q);
        plot(r_extrap, q_extrap);
        title('Evolution of q.\nabla \theta fits over time');
        xlabel('r'); ylabel('q.\nabla \theta');
        %}
        %pause;
        
        
        %Find predicted q(a)
        a_i = round(a_n(1)/dr_extrap) + 1;
        q_a =  q_extrap(1, a_i) * ones(z_num, 1);
        
    end
    
    %Enforce relaxation
    lambda = constants('relaxation_lambda');
    %a_new = a_n + dt*lambda*(a_extrap - a_n);
    %a_new = a_n + dt*lambda*(a_interp - a_n);
    a_new = a_n + dt * lambda * q_a;
      
    a_new = min(max(a_new, 0), b);
 
    
end

fprintf('new a: %1.6f, and b = %1.6f \n', a_new(1), b);

end
