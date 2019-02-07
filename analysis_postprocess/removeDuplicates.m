function removeDuplicates(Rms, Rs, Hs)

points = 40;

for Rm_i = 1:numel(Rms)
    Rm = Rms(Rm_i);
    for R_i = 1:numel(Rs)
        R =Rs(R_i);
        for H_i = 1:numel(Hs)
            H = Hs(H_i);

            points = 40;
            %folder = [pwd '/data/'];
            folder = [pwd '/data/theta_inf_1.4/'];

            %file = strcat('mushyLayerPrevSteadyState_relaxed_Points',num2str(points),'Rm',num2str(Rm),'H',num2str(H)','R',num2str(R),'*.mat');
            file = ['mushyLayerPrevSteadyState_relaxed_Points' num2str(points) 'Rm' num2str(Rm) 'H' num2str(H) 'R' num2str(R) '*.mat'];
            filesearch = [folder file];
            filelist = dir(filesearch);

            if numel(filelist) > 1

                %1. Find the newest file
                i_newest = 1;
                datenum_newest = filelist(1).datenum;

                for i = 2:numel(filelist)
                    if filelist(i).datenum > datenum_newest;
                        datenum_newest = filelist(i).datenum;  
                        i_newest = i;
                    end
                end

                %Delete all but the newest
                for i = 1:numel(filelist)
                    if filelist(i).datenum < datenum_newest;
                        delete([folder filelist(i).name]);
                        fprintf('Deleted %s\n', [folder filelist(i).name]);
                    end
                end
                fprintf('Kept %s\n', [folder filelist(i_newest).name]);
            end
            
        end
    end
end
