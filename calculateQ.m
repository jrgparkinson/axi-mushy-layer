function q = calculateQ(r,z,theta,psi)

[dr, dz, ~,~] = meshGridProperties(r,z);

[theta_r, theta_z] = gradient2order(theta, dr, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r, z, psi);
q_r = u_r;  q_z = u_z-1;
q = theta_r .* q_r + theta_z.*q_z;