% Code to Partition CSV files into 2D emission correlation slices
clc
clear
tic

main_folder = pwd;
addpath(main_folder);

%% 
sourceflag = 'raw_data'
% "monte/data"

x = input( ['you have selected data source as : ' sourceflag ' press enter to continue or press Q to quit]']);
if x=='Q'
    error('partitioning exit manually')
end

switch sourceflag
    case 'monte'
       [file_name,data_path] = uigetfile; 
       load_dat = load([data_path file_name ]);
       disp 'data loaded'
        
    case 'raw_data'
        file_name = 'mau';
        data_path = uigetdir; 
        load_dat = 0;
    otherwise
        error('choose "monte" or "data"')
end

[~,name_str] = fileparts(data_path);
%% Partition windows

p1 = [1E-6 10E-6];
p2 = [10E-6 20E-6];
p3 = [20E-6 30E-6];
p4 = [30E-6 40E-6];
p5 = [40E-6 50E-6];
p6 = [50E-6 60E-6];
p7 = [60E-6 70E-6];
p8 = [70E-6 80E-6];
p9 = [80E-6 90E-6];
p10 = [90E-6 100E-6];
p11 = [100E-6 200E-6];
p12 = [200E-6 300E-6];
p13 = [300E-6 400E-6];
p14 = [400E-6 500E-6];
p15 = [500E-6 600E-6];
p16 = [600E-6 700E-6];
p17 = [700E-6 800E-6];
p18 = [800E-6 900E-6];
p19 = [900E-6 1E-3];
p20 = [1E-3 2E-3];
p21 = [2E-3 5E-3];
p22 = [5E-3 10E-3];
p23 = [10E-3 20E-3];
p24 = [20E-3 50E-3];
p25 = [50E-3 100E-3];
p26 = [100E-3 200E-3];
p27 = [200E-3 300E-3];
p28 = [300E-3 500E-3];
p29 = [500E-3 700E-3];
p30 = [700E-3 1];
p31 = [1 2];
p32 = [2 3];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% frames = p31;
frames = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20;p21;p22;p23;p24;p25;p26;p27;p28;p29;p30;p31;p32];
%     frames = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15];
%      frames = [p33; p34; p35; p36 ;p37; p38; p39; p40]

%% parfor loop
fid =(11:length(frames)+10);
% fid = 11

parfor i = 1:length(fid)
    frames(i,:)
    loop_window_lb = frames(i,1);
    loop_window_ub = frames(i,2);
    fname_str = [name_str '-' num2str(fid(i))];
   main_exp{i}.partition =  all_frames_partition(fname_str,loop_window_lb,loop_window_ub,sourceflag,data_path,file_name,load_dat,main_folder)
   cd(main_folder)
end
toc

data_folder = GetParentFolder(data_path);
[~,part_str] = fileparts(data_folder);
part_dir = 'Partitioned_data';
cd(data_folder)
mkdir(part_dir)
save_file_name = ['partitioned_' part_str];
save(save_file_name , 'main_exp','-v7.3') 
movefile([save_file_name '.mat'], part_dir);
cd(main_folder)

fprintf('*****************partitioning complete ***************** \n')
