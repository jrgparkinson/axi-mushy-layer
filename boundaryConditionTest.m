%Test that the steady state satisfies the boundary conditions

%Test for different grid spacings
r_nums = [25 30 35 40 45 50];
for i = 1:numel(r_nums);
    r_num = r_nums(i);
    filename = strcat('mushyLayerPrevSteadyState',num2str(r_num),'.mat');
    
    %Load r, z, psi, theta, a, C, constants
    load(filename);
    Da = constants('Da');
    Rm = constants('Rm');
    b = constants('b');
    
    [dr, dz, r_num, z_num] = meshGridProperties(r,z);
    %{
    psi = psi;
    plot(z(:, 1), psi(: , 1));
    xlabel('r'); title('\psi at r=b'); ylabel('\psi');
    pause;
    plot(z(:, 1), theta(: , 1));
    xlabel('r'); title('\theta at r=b'); ylabel('\theta');
    pause;
    %}
    
    %Fixed boundary conditions
    psi_right = psi(:, end); %Should = 0
    
    psi_bottom = psi(1, :); %Should = 0
    theta_bottom = theta(1, :); %Should = -1
    
    theta_top = theta(end, :); %Should = 1
    
    %Simple dirichlet boundary conditions
    psi_z_top = (-3*psi(end, :) + 4*psi(end-1, :) - psi(end-2, :))/(2*dz); %Should be 0
    
    theta_r_right = (-3*theta(:, end) + 4*theta(:, end-1) - theta(:, end-2))/(2*dr); %Should be 0
    
    %Chimney
    %psi_r_left = (-3*psi(:, 1) + 4*psi(:, 2) - psi(:, 3))/(2*dr);
    %theta_r_left = (-3*theta(:, 1) + 4*theta(:, 2) - theta(:, 3))/(2*dr);
    [psi_r, psi_z] = gradient2order(psi, dr, dz);
    [theta_r, theta_z] = gradient2order(theta, dr, dz);
    
    heat_test = b*theta_r(:, 1) - ...
        (psi(:, 1)-0.5*b^2) .* theta_z(:, 1);  %Should be zero
    
    psi_test = abs(psi(:, 1) - ...
        (a.^4)/(16*Da) .* (psi_r(:,1)/b + Rm*(theta(:, 1) - C(:, 1))) + 0.5*b*psi_r(:, 1)); %Should be zero
    
    %Add 0.1 do avoid dividing by very small things
    heat_test_frac = heat_test./(b*theta_r(:, 1) + 1);
    psi_test_frac = psi_test./(psi(:, 1) + 1);
    
    plot(z, heat_test_frac); xlabel('z');  title('fractional error in heat bc at chimney');
    pause; 
    plot(z, psi_test_frac); xlabel('z');  title('fractional error in \psi bc at chimney');
    pause;
    
    %Make metrics
    heat_max_metric(i) = max(heat_test_frac);
    heat_av_metric(i) = nanmean(heat_test_frac);
    
    psi_max_metric(i) = max(psi_test_frac);
    psi_av_metric(i) = nanmean(psi_test_frac);
    
end

plot(r_nums.^2, heat_max_metric);
xlabel('(\Delta r)^2');
title('Max fractional error in heat bc at chimney');
pause;
plot(r_nums.^2, heat_av_metric);
xlabel('(\Delta r)^2');
title('Average fractional error in heat bc at chimney');
pause

plot(r_nums.^-2, psi_max_metric);
xlabel('(\Delta r)^2');
title('Max fractional error in \psi bc at chimney');
pause;
plot(r_nums.^-2, psi_av_metric);
xlabel('(\Delta r)^2');
title('Average fractional error in \psi bc at chimney');
pause


