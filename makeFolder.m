function [save_dir] = makeFolder(mat_fileName, base_name, folder_name)

% 1. Make new folder
unix_time = num2str(round(posixtime(datetime)));
save_dir = [folder_name '\' base_name '_' unix_time];
system(['mkdir ' save_dir]);
disp(save_dir); % Display title of directory
save_dir = [save_dir '\'];

% 2. Save copy of running script to new folder
mat_newName = [save_dir '\' mat_fileName unix_time '.m'];
cpstring = [mat_fileName '.m ' mat_newName];
[~,~] = system(['copy ', cpstring]); % Copy [matlab file] to [new matlab file in new folder]

end
