function [dr, dz, r_num, z_num] = meshGridProperties(r, z)

lengths = size(r);
z_num = lengths(1);
r_num = lengths(2);

dr = (r(1, end) - r(1, 1) ) / (r_num - 1);
dz = (z(end, 1) - z(1, 1) ) / (z_num - 1);


end