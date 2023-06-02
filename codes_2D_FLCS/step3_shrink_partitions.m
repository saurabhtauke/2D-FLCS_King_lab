% code to shrink the partions into 256 x 256 sizes
% There is an opportunity to integrate this step into the partitioning
% algorithm itself

clc
clear
main_folder = pwd;
addpath(main_folder);   

%% load data
[select_file, file_path] = uigetfile('*.mat');
load([file_path select_file]);

num_parts  = max(size(main_exp));
ndata = max(size(main_exp{1}.partition.TwoD_d1d1));

% specify shrink parameters
shrink_to = 256;
shrink_by = ndata/shrink_to;

new_part = cell(size(main_exp));

% loop over all partitions 
for i = 1:num_parts
    
    OneD_dec1 = main_exp{i}.partition.OneD_dec1;
    sh1D_dec1 = shrink1D(OneD_dec1,shrink_by);

    OneD_dec2 = main_exp{i}.partition.OneD_dec2;
    sh1D_dec2 = shrink1D(OneD_dec2,shrink_by);    
    
    tot2D = main_exp{i}.partition.TwoD_d1d1+main_exp{i}.partition.TwoD_d1d2+main_exp{i}.partition.TwoD_d2d1+main_exp{i}.partition.TwoD_d2d2;
    shData2D = shrink2D(tot2D,shrink_by,ndata);
    
    partition.sh1D_dec1 = sh1D_dec1;
    partition.sh1D_dec2 = sh1D_dec2;
    partition.sh2D = shData2D;
    new_part{i} = partition;
    
end

% store new partitions
cd(file_path)
name_str = ['shrinked_new_' select_file];
save(name_str, 'new_part')
cd(main_folder)

fprintf('*****shrinking partitions complete ******** \n')

%% functions

function [sh_trace] = shrink1D(in_trace,shrink_by)
    % data binning by a specified factor
    
    dat_size = length(in_trace);
    shrinked_size = ceil(dat_size/shrink_by);
    shrink = 1:shrink_by:dat_size;
    
    sh_trace = zeros(shrinked_size,1);

    for i =1:shrinked_size-1
        sh_trace(i) = sum(in_trace(shrink(i):(shrink(i+1)-1)));

    end
    sh_trace(shrinked_size) = sum(in_trace(shrink(shrinked_size):length(in_trace)));
    
end


function [shData2D] = shrink2D(in_2D_data,shrink_by,dat_size)
    
    shrinked_size = ceil(dat_size/shrink_by);  
    shrink = 1:shrink_by:dat_size;
    shrink_datvec = zeros(shrinked_size,shrinked_size);
    
    for i = 1:shrinked_size-1
        for j = 1:shrinked_size-1
            shrink_datvec(i,j) = sum(sum(in_2D_data(shrink(i):(shrink(i+1)-1),shrink(j):(shrink(j+1)-1))));
        end
    end      
    
    shData2D = shrink_datvec;  

end
