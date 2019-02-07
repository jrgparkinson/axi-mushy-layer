function theta_inf = calculateThetaInfinity(r, z, theta, psi)

[dr, dz, ~, ~] = meshGridProperties(r,z);

[~, theta_z] = gradient2order(theta, dr, dz);
[psi_r, ~] = gradient2order(psi, dr, dz);


theta_inf = theta_z(end, end) / (1-psi_r(end, end)./r(end, end));