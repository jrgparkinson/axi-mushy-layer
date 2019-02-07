%load('Rm40H0.5R0.4a0.04b0.05Interp.mat');
load('Rm40H0.5R0.4a0.04b0.05.mat');
%load('Rm40H0.5R0.4a0.04b0.05_oldC.mat');

plot(z(:, 1), C(:, 1)); pause;

[dr, dz, r_num, z_num] = meshGridProperties(r,z);
j = round(0.5*z_num);

[theta_r, theta_z] = gradient2order(theta, dr, dz);
[u_r, u_z] = axiVelocitiesFromPsi(r, z, psi);

q = u_r.*theta_r + (u_z-1).*theta_z;

clf;
figure(1); hold on;
plot(r(j, :), q(j, :));
pause;

q_4 = q(j, 1:4).';
r_4 = r(j, 1:4).';

% Fit to a cubic curve
fo = fitoptions('Method','NonlinearLeastSquares',...
    'StartPoint',[-100 10 0.006 0.05]);

%f = fittype('(a*x^3+b*x^2+c*x+d)/x','options', fo);
f = fittype('a*x^2+b*x+c+d/x','options', fo);
fitobject = fit(r_4, q, f);
coeffs = coeffvalues(fitobject);

plot(fitobject, r_4, q);

a1 = coeffs(1); a2=coeffs(2); a3=coeffs(3); a4 = coeffs(4);

r_extrap = 0:0.001:b;
%q_extrap = a1*r_extrap.^2+a2*r_extrap+a3+a4./r_extrap;
q_extrap = (a1*r_extrap.^3+a2*r_extrap.^2+a3*r_extrap+a4)./r_extrap;
plot(r_extrap, q_extrap);
pause;