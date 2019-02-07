function Flux = calculateSolidHeatFlux(r,z,theta,psi, constants,a, C, scaled)

[dr, dz, r_num, z_num] = meshGridProperties(r,z);

[theta_r, theta_z] = gradient2order(theta, dr, dz);
[psi_r, psi_z] = gradient2order(psi(:,:), dr, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r, z, psi);
q_z = u_z - 1;
q_r = u_r;

heat_r = q_r .* theta - theta_r; %heat in transport in r
heat_z = q_z .* theta - theta_z; %heat transport in z

%[r_tot, z_tot, theta_tot, psi_tot, Q, i_min, b_i] = extrapolateSolution(r, z, theta, psi, a, C, constants);
%[dr_tot, dz_tot, ~,~] = meshGridProperties(r_tot,z_tot);
%[psi_tot_r, psi_tot_z] = gradient2order(psi_tot, dr_tot, dz_tot);

b = constants('b');
R = constants('R');
H = constants('H');

z_ms = 1; 
z_mo = z_num;
r_b = 1; %r index for r=b

z_exclude_singularity = round(z_mo*0.9);
%z_exclude_singularity = z_mo;

r_exclude_singularity = round(r_num*0.1);
%r_exclude_singularity = r_b;

%Solid flux
integral_ms = trapz(r(z_ms, r_b:end), heat_z(z_ms, r_b:end));
F_ms = integral_ms;

%Chimney flux (approx - calculate at b not a)
integral_mc = trapz(z(1:z_exclude_singularity, r_b), heat_r(1:z_exclude_singularity, r_b));
F_mc = integral_mc; %(2*pi*H)

%Ocean flux
integral_mo = trapz(r(z_mo, r_exclude_singularity:end), heat_z(z_mo, r_exclude_singularity:end));
F_mo = integral_mo;


%F_chimney_centre = trapz(z_tot(:, 1), psi_tot_z(:, 1).*(theta_tot(:, 1) - 1));

fprintf('F_mc/R: %1.3f, F_mo/R: %1.3f, F_ms/R: %1.3f \n', F_mc/R, F_mo/R, F_ms/R);

Flux = F_ms;


if scaled == 1
    Flux = Flux/(R);
end

Flux = abs(Flux/pi);