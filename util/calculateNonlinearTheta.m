function variation = calculateNonlinearTheta(r, z, H, theta, r_range)

[dr, dz, r_num, z_num] = meshGridProperties(r,z);

variation = nan*r_range;

theta_linear = -1 + z/H;
integrand = abs(theta - theta_linear);

frac_error = abs(integrand./(theta -0.01));


for r_i = r_range
    %variation(r_i) = trapz(z(:, r_i), integrand(:, r_i));
    %variation(r_i) = mean(integrand(:, r_i));
    
    variation(r_i) = mean(frac_error(:, r_i));
end

%variation = trapz(z(:, r_num), integrand(:, r_num));

%average_deviation = mean(integrand(:, r_num));
%fprintf('%1.5f \n', variation(r_num));

%mean_av = average_deviation * ones(z_num);

%plot(z(:, r_num), integrand(:, r_num)./ (z(:, r_i) + 0.1), z(:, r_num), mean_av); pause;

%plot(z(:, r_num), theta_linear(:, r_num),z(:, r_num),theta(:, r_num) ); pause;
%plot(z(:, r_num), abs(integrand(:, r_num)./(theta(:, r_num) -0.01))); pause;

end