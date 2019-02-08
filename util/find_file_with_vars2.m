function filename = find_file_with_vars(R,H,Rm,points, flst)

filename_search = strcat('mushyLayerPrevSteadyState_relaxed_Points',num2str(points),'Rm',num2str(Rm),'H',num2str(H),'R',num2str(R),'b');

ix=regexp(flst,filename_search);
ix=~cellfun('isempty',ix);
matched_files=flst(ix);

if numel(matched_files) == 0
    filename = false;
elseif numel(matched_files) == 1
    filename = strcat('/data/', matched_files{1});
else
    filename = strcat('/data/', matched_files{1});
    fprintf('Error - found multiple file names for the search: \n \t  %s \n', filename_search);
end



end