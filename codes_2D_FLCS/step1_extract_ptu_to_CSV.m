% Code to 
% 1. extract ptu files
% 2. convert ptu to csv

clc
clear 
main_folder = pwd;
addpath(main_folder);

%% load folder

% select folder containing data
folder_path = uigetdir();
cd(folder_path);
pathname = [folder_path '\'];

% list all .ptu files
ptu_files = ls('*ptu');

%% run ptu_process in loop

loop_len = length(ptu_files(:,1));

for i = 1:loop_len
    filename = ptu_files(i,:);     
    ptu_read_in_loop(filename,pathname);
    wt =  waitbar(i/loop_len);
end 
close(wt)

%% save to CSV folder

csv_dir = 'CSV_converted';
mkdir(csv_dir)
csv_list = ls('*csv');  
for i  = 1:length(csv_list(:,1))
    movefile(csv_list(i,:), csv_dir);
end

cd(main_folder);
fprintf('Conversion  completed \n')