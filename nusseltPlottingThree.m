close all;
clear all;

Ras = [10 4*pi^2 50 100 200 250 300];
Nu_analytic = [1  1  1.45 2.651 3.813 4.199 4.523 ];
Nu_calculated_numeric = [1 1 1.45928 2.69869 3.87044 4.27359 4.61144];
Nu_calculated_matlab = [1 1 1.44844 2.63780 3.80152 4.18411 4.50164];

psi_analytic = [0 0 2.112 5.377 8.942 10.253 11.405];
psi_calculated = [0 0 2.1086 5.36441 8.93465 10.24223  11.38452];

plot(Ras,  Nu_analytic,'or', Ras, Nu_calculated_matlab, '--+b');
xlabel('Ra');
ylabel('Nu');
axis([0 350 0 5]);
pause;
%{
plot(Ras,  psi_analytic,'--*g', Ras, psi_calculated, '--*b');
%title('Maximum psi calculated by Caltagirone (green) and me (blue)');
xlabel('Ra');
ylabel('\psi');
axis([0 350 0 13]);
pause;
%}
%{
%Difference between solutions
Nu_frac_diff = (Nu_calculated - Nu_analytic) ./ Nu_analytic;
plot(Ras,  Nu_frac_diff, '--*');
title('Nusselt number fractional difference');
xlabel('Fractional difference');
ylabel('Nu');
pause;


%Convergence for Ra = 200
dxs = [1/32 1/48 1/50 1/60 1/70];
Nu = [3.68123 3.72999 3.73361 3.74765 3.75729];
Nu_frac_diff = abs(Nu - 3.813)/3.813;
plot(dxs, Nu_frac_diff);
title('Fractional difference in Nu(Ra=200) for \epsilon = 0.1, \Delta t = 0.01');
xlabel('\Delta x');
ylabel('(\Delta Nu)/(Nu_{Caltagirone})');
pause;
%}
