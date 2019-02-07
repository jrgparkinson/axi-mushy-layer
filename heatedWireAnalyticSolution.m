% Calculate analytic solution
function [T, psi, psi_z, psi_r] = heatedWireAnalyticSolution(r_points, z_points, a_nc, T_inf, Rm, options, alpha_m)

n_max = 15;

[f, theta, n_points] = solveODE_bvp(a_nc, n_max, options);
%[f_ode, theta_ode, n_points_ode] = solveODE(a_nc, n_max);

%plot(n_points, f_arr);

%r_num = numel(r_points);        z_num = numel(z_points);
[dr, dz, r_num, z_num] = meshGridProperties(r_points, z_points);

%T = zeros(r_num, z_num);   psi = zeros(r_num, z_num);
T = 0*meshgrid(1:r_num, 1:z_num);     psi = 0*meshgrid(1:r_num, 1:z_num);
psi_r = 0*meshgrid(1:r_num, 1:z_num); psi_z = 0*meshgrid(1:r_num, 1:z_num);

for r_i = 1:r_num
    r = r_points(1, r_i);
    n = Rm * r^2;
    
    %z_j = 1 boundary
    %T(r_i, 1) = T_inf;
    %psi(r_i, 1) = 0;
    %psi_r(r_i, 1) = 0;
    
    T(1, r_i) = T_inf;
    psi(1, r_i) = 0;
    psi_r(1, r_i) = 0;
    
    for z_j = 2:z_num
        z = z_points(z_j, 1);
        T_w = T_inf + z;
        
        if (n < a_nc)
            %T(r_i, z_j) = T_w;
            %psi(r_i, z_j) = 0;
            %psi_z(r_i, z_j) = 0;
            %psi_r(r_i, z_j) = 0;
            
            T(z_j, r_i) = T_w;
            psi(z_j, r_i) = 0;
            psi_z(z_j, r_i) = 0;
            psi_r(z_j, r_i) = 0;
            
        elseif n > n_max
            T(z_j, r_i) = T_inf;
            psi(z_j, r_i) = psi(z_j, r_i-1);
            
            psi_r(z_j, r_i) = 0;
            psi_z(z_j, r_i) = psi_z(z_j, r_i - 1);
            
        else
            f_n = interp1(n_points, f, n);
            theta_n = interp1(n_points, theta, n);
            
            psi(z_j, r_i) = alpha_m * z * f_n;
            T(z_j, r_i) = T_inf + z*theta_n;
            psi_z(z_j, r_i) = f_n;
            psi_r(z_j, r_i) = Rm * z * r * theta_n;
        end
        
    end
    
end


end

function fp=G(n, f)
fp=zeros(3,1);
fp(1) = (-f(1) + 0.5*f(3)*f(1) + 0.5*f(3)*f(2));
fp(2) = f(1);
fp(3) = f(2);
end

function [f, theta, n_points] = solveODE_bvp(a_nc, n_max, options)

solinit = bvpinit(linspace(a_nc,n_max,200), @guess);
sol = bvp4c(@func,@bcs,solinit,options);

% Plot t v.s. y1(t):figure;
f = sol.y(1,:);
theta = 2*sol.y(2,:);
n_points = sol.x;
%figure; plot( sol.x, theta, '-bo' ); grid on;
%title('\theta');
end

function [v] = guess(x)

v = [ x;
    exp(-x);
    -exp(-x)];
end

function [res] = bcs(ya,yb,beta)
% BCS -
%
res = [ ya(1) - 1;
    ya(2) - 0.5;
    yb(2) - 0;
    ];
end

function [ ode ] = func(x,y)
% F -
% y(1)=F(x); y(2)=F'(x); y(3)=G(x);
%-2*(x*y(1)-y(3))*y(2) - 2*(1-y(1)^2);
%x*y(2);

ode = [
    y(2);
    y(3);
    (y(2)^2 - y(1)*y(3) - 2*y(3))/(2*x);
    ];
end