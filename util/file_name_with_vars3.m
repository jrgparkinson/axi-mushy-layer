function filelocation = file_name_with_vars3(R, Rm, theta_inf)

points = 40;

%Get all files with the right theta infinity, R and Rm

folder = [pwd '/data/theta_inf_' num2str(theta_inf) '/'];

file = strcat('mushyLayerPrevSteadyState_relaxed_Points',num2str(points),'Rm',num2str(Rm),'H*R',num2str(R),'b*.mat'); 
filesearch = [folder file];
filelist = dir(filesearch);

if numel(filelist) > 0
    
    if numel(filelist) > 1
        %Find the file with the closest theta_inf to 1.4 (and delete the
        %others)
        theta_diff_smallest = 0.02;
        i_smallest = 1;
        
        for i = 1:numel(filelist)
            load([folder filelist(1).name]);
            theta_infinity_file = calculateThetaInfinity(r,z,theta,psi(:,:));
            theta_diff(i) = abs(theta_infinity_file - theta_inf)/theta_inf;
            
            if theta_diff(i) < theta_diff_smallest
                theta_diff_smallest = theta_diff(i);
                i_smallest = i;
            end
        end
        
        for i = 1:numel(filelist)
            if theta_diff(i) > theta_diff_smallest 
                delete([folder filelist(1).name]);
            end
        end
        
        filelocation = [folder filelist(i_smallest).name];
        
    else
        %Only one file found
        filelocation = [folder filelist(1).name];
    end
    
else
    filelocation = false;
end
