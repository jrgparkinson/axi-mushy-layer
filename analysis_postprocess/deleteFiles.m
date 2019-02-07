function deleteFiles

R_vals = 0.268:0.003:0.28;%0.01:0.003:0.45; 
Rm_vals = [65]; %45:0.5:80; 

R_num = numel(R_vals);  Rm_num = numel(Rm_vals);

%flst=dir([pwd '/data/*.mat']);
%files_to_search={flst.name};

%delete('./data/mushyLayerPrevSteadyState_relaxed_Points40Rm60H0.469R0.316b0.050525.mat');

for R_i = 1:R_num
    R = R_vals(R_i);
        for Rm_i = 1:Rm_num
            Rm = Rm_vals(Rm_i);
            
            
            filename = strcat('mushyLayerPrevSteadyState_relaxed_Points40Rm',num2str(Rm),'H*R',num2str(R),'b*.mat');
            filelocation = strcat('./data/', filename);
            theta_inf_location = strcat('./data/theta_inf_1.4/', filename);
            %filename = './data/mushyLayerPrevSteadyState_relaxed_Points40Rm60H*R0.316b*.mat';
            delete(filelocation);
            delete(theta_inf_location);
            fprintf('Deleted %s \n', filename);
            
        end

end