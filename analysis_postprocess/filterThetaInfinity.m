%Get all files and calculate theta infinity for them
function filterThetaInfinity(theta_inf)

flst=dir([pwd '/data/mushyLayerPrevSteadyState_relaxed_Points40*.mat']);
files_to_search={flst.name};

for i=1:numel(files_to_search)
    filename =  files_to_search{i};
    filelocation = [pwd '/data/' filename];
    
    load(filelocation);
    theta_inf_file = calculateThetaInfinity(r,z,theta,psi(:,:));
    theta_diff = (theta_inf_file - theta_inf)/theta_inf;
    
    if abs(theta_diff) < 0.02
        copyto =strcat('./data/theta_inf_',num2str(theta_inf),'/',filename);
        copyfrom = filelocation;
        copyfile(copyfrom, copyto,'f');
    end
    
    fprintf('Parsed %d/%d files \n', i, numel(files_to_search));
end

