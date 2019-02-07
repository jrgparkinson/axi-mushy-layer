vars = 'Rm40H0.5R0.4a0.04b0.05';
load(strcat('steadyState',vars,'.mat'));


figure;
mesh(r, z, theta); 
xlabel('r'); ylabel('z'); title('\theta steady state extrapolated to r=0');
savefig(strcat('thetaSteadyState',vars,'.fig'));

figure;
mesh(r, z, psi); 
xlabel('r'); ylabel('z'); title('\psi steady state extrapolated to r=0');
savefig(strcat('psiSteadyState',vars,'.fig'));

figure;
mesh(r, z, Q_m); 
xlabel('r'); ylabel('z'); title('Q = q.\nabla(\theta) steady state extrapolated to r=0');
savefig(strcat('QSteadyState',vars,'.fig'));

