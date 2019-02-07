%Extend solution to a < r < b

function [r_tot, z_tot, T_tot, psi_tot, Q, i_min, b_i] = extrapolateSolution(r, z, T, psi, a_final, C, constants)

b = constants('b');
Rm = constants('Rm');
Da = constants('Da');
R = constants('R');
H = constants('H');

[dr, dz, r_num, z_num] = meshGridProperties(r,z);
z_min = z(1,1);
z_max = z(end, 1);

%Need many grid points in r direction to resolve chimney
dr_extrap = b/50;

%Sample entire grid at dr_extrap
[r_new, z_new] = meshgrid(b:dr_extrap:R, z_min:dz:z_max);
T_new = interp2(r, z, T, r_new, z_new, 'spline');
psi_new = interp2(r,z,psi, r_new, z_new, 'spline');

[r_extrap, z_extrap] = meshgrid(0:dr_extrap:(b-dr_extrap), z_min:dz:z_max);
r_tot = cat(2, r_extrap, r_new);    z_tot = cat(2, z_extrap, z_new);
T_extrap = nan * r_extrap;  psi_extrap = nan * r_extrap;
T_tot = cat(2, T_extrap, T_new);    psi_tot = cat(2, psi_extrap, psi_new);
T_m = T_tot;    psi_m = psi_tot;

C_tot = interp1(z(:, 1), C, z_tot(:, 1));

s = size(r_extrap); b_i = s(2)+1;  %b_i is the index of r=b

theta_r_b = (- 3 * T_tot(:, b_i) + 4 * T_tot(:, b_i+1) - T_tot(:, b_i+2)) / (2*dr_extrap);
psi_r_b = (- 3 * psi_tot(:, b_i) + 4 * psi_tot(:, b_i+1) - psi_tot(:, b_i+2)) / (2*dr_extrap);

% 0th order approx for psi
%B = a_final(:, 1).^4/(16*Da) .* (psi_r_b/b + Rm*(T_tot(:, b_i) - C(:, 1)) );

% 1st order approx for psi
%B = a_final(:, 1).^4/(16*Da) .* (psi_r_b/b + Rm*(T_tot(:, b_i) - C(:, 1)) ) + (1/6) * (b^3 - a_final.^3) * Rm .* theta_r_b(:, 1);

a = a_final(1);

for j = 1:z_num
    
    a_i = round(a/dr_extrap) + 1; % index of the point r=a
    %i_min = round(a/dr_extrap);
    i_min = 1;

    for i = (b_i-1):-1:i_min 
        r_i = r_tot(1, i);
        r_iPlusOne = r_tot(1, i+1);
        if (i >= a_i)
            % 0th order approximation for T
            %T_concat(j, i) = T_final(j, 1);

            % 1st order approx for T
            T_tot(j, i) = T_tot(j, b_i) - (b-r_tot(j, i)) * theta_r_b(j, 1);
            
        else
            T_tot(j, i) = T_tot(j, a_i);
        end  
        
         
        X = Rm*b*theta_r_b(j) + psi_r_b(j)/b;
        
        T_m(j, i) = T_m(j, b_i) - (b-r_tot(j, i)) * theta_r_b(j, 1);
        psi_m(j,i) = psi_tot(j, b_i) + 0.5 * (r_i^2 - b^2) * X - (1/3)*(r_i^3 - b^3)*Rm*theta_r_b(j);
        
        
        % approx for psi
        %psi_tot(j, i) = psi_tot(j, i+2) + 4*dr_extrap/(r_iPlusOne) * (-psi_tot(j, i+1) + B(j));
        
        if i >= a_i 
  
            psi_tot(j, i) = psi_tot(j, b_i) + 0.5 * (r_i^2 - b^2) * X - (1/3)*(r_i^3 - b^3)*Rm*theta_r_b(j);
            
        else 
            psi_r_a = (- 3 * psi_tot(j, a_i) + 4 * psi_tot(j, a_i+1) - psi_tot(j, a_i+2)) / (2*dr_extrap);
            p_z = - psi_r_a/(a*Rm) - T_tot(j, a_i);
            Const = (Rm/Da) * (p_z + C_tot(j, 1));
            
            A = (1/16) * Const;
            B = (1/2) * (psi_r_a/a - 0.25 * Const*a^2);

            %approx = A * (r_i^4-a^4) + 0.5*B*(r_i^2-a^4);
            %psi_tot_a = psi_tot(j, a_i);
            
            %initial_diff =  psi_tot(j, a_i) - (A * (a^4) + B*(a^2));
            
            %psi_tot(j, i) = psi_tot(j, a_i) + A * (r_i^4-a^4) + B*(r_i^2-a^2);
            psi_tot(j, i) = A * (r_i^4) + B*(r_i^2);
            
            
        end
        %psi_tot(j, i+2) + 4*dr_extrap/(r_iPlusOne) * (-psi_tot(j, i+1) + B(j));
        
    end
end
        
%Calculate derivatives
%{
[T_r, T_z] = gradient2order(T_tot, dr_extrap, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r_tot, z_tot, psi_tot);

q_r = u_r;      q_z = u_z-1;
Q = q_r.*T_r + q_z.*T_z;
%}
Q = calculateQ(r_tot, z_tot, T_tot, psi_tot);

%{
[T_r_m, T_z_m] = gradient2order(T_m, dr_extrap, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r_tot, z_tot, psi_m);
q_r = u_r;  q_z = u_z-1;
%}
Q_m = calculateQ(r_tot, z_tot, T_m, psi_m);

theta = T_m; psi = psi_m;   r = r_tot;  z = z_tot;  

vars = strcat('Rm',num2str(Rm),'H',num2str(H),'R',num2str(R),'a',num2str(a),'b',num2str(b));
title_vars = strcat('Rm=',num2str(Rm),', H=',num2str(H),', R=',num2str(R),', a=',num2str(a),', b=',num2str(b));

%{
filename = strcat(vars,'.mat');
a = a_final;
save(filename, 'r', 'z', 'theta', 'psi', 'Q_m', 'C', 'a');
fprintf('Saved %s \n', filename);
%}

end