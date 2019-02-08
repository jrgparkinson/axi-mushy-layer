function H = find_H(theta_inf, R, Rm, r_num, files_to_search)
H_test = 0.2:1e-3:0.3; %1e-5 steps
H_num = numel(H_test);

for H_i = 1:H_num
    H = H_test(H_i);
    filename = find_file_with_vars(R,H,Rm, r_num, files_to_search);
    
    if filename ~= false
        
        load([pwd filename]);
        theta_infinity = calculateThetaInfinity(r,z,theta,psi(:,:));
        theta_diff = (theta_infinity - theta_inf)/theta_inf;

        if abs(theta_diff) < 0.01
            %We've found H = return it
            fprintf('Rm: %2.2f, R:%1.5f, H: %1.5f, theta inf: %1.5f, theta_diff: %1.5f ** \n', Rm, R, H, theta_infinity, theta_diff);
            break;
        elseif abs(theta_diff) < 0.05
            fprintf('Rm: %2.2f, R:%1.5f, H: %1.5f, theta inf: %1.5f, theta_diff: %1.5f -- \n', Rm, R, H, theta_infinity, theta_diff);
        elseif abs(theta_diff) < 0.3
            fprintf('Rm: %2.2f, R:%1.5f, H: %1.5f, theta inf: %1.5f, theta_diff: %1.5f - \n', Rm, R, H, theta_infinity, theta_diff);            
        end
    end
    
    % If we didn't find H, return nan
    H = nan;
    theta_infinity=nan;
end