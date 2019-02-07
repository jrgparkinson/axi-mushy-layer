function [u_r, u_z] = axiVelocitiesFromPsi(r, z, psi)
[dr, dz, ~, ~] = meshGridProperties(r, z);

[psi_r, psi_z] = gradient2order(psi, dr, dz);

u_r = - psi_z./r;
u_z = psi_r./r;

end
