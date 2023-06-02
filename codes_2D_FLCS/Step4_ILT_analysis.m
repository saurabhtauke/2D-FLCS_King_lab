%% 2D ILT code
% Check optimum IRF shift by specifying a range of integers in
% IRF_1D_shift_range

%%
clc
clearvars -except IJM main_exp select_file  % Clear all variables
tic
main_folder = 'E:\Saurabh\TUTORIAL_2DFLCS\codes_2D_FLCS';
addpath(main_folder);

%% specify inputs

str_name = 'sample';

IRF_1D_shift_range = -267;

InitTime_T1 = 0.1;
FinalTime_T1 = 10;                         
L_basis = 50;


%%
new_dir =  [str_name '-batch_processed']
if isfile(['analyzed_' str_name '.mat'])
%     error( ' ERROR : file exists; delete the existing files or rename the analysis file')
      warning( ' ERROR : file exists; delete the existing files or rename the analysis file')
end
   
   fit_1d = 1;
   fit_2d = 1;

if length(IRF_1D_shift_range)>1
    fit_2d = 0;
end
IRF_2D_shift = IRF_1D_shift_range;

%% global fit_parameters

global microtime_resolution
global ndata

microtime_resolution = 0.004;
shrink_by = 16;

Singular_vals = 44;
cond_num_bound = 1E9; % 1E6 is ok

input_alpha = 0;
   
   ndata = 4093;
%    error('FIX: implicit ndata determination')

   IRF_flag = 1;
  
   fixed_basis_flag = 0;
   if  fixed_basis_flag
       L_basis = 3;
   end
   
   nrepeat = 1; 
   noise_amplitude1d = 0E-3;
   noise_amplitude2d = 0E-4;
   
    basis_scale = 'log';
    
%% load IRF
[nIRF,IRF_time]= load_IRF('irf_wavelength_p3_40mhz.mat');
      
%% loop data folder

if exist('new_part','var')
    %chill
    warning('using pre loaded main_exp')
    tmp = input('press any key to continue  ','s');
else
    [select_file, file_path] = uigetfile('*.mat');
    load([file_path select_file])
end

dat_fold = GetParentFolder(GetParentFolder(file_path));
cd(dat_fold)

if isfolder(new_dir)
%     error('ERROR:folder already exists: RENAME NEW_DIR')
else
    mkdir(new_dir);
end

n_files = length(new_part);

ILT_frame = zeros(L_basis,L_basis,n_files);
pd = zeros(L_basis,L_basis,n_files);

%% analysis loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make loop from here
if fit_1d
for MAIN_LOOP = 1:1

cd(main_folder)

data1D = new_part{MAIN_LOOP}.sh1D_dec1;
data1D = process_data(data1D);

% data2D_all = main_exp{MAIN_LOOP}.partition.TwoD_d1d1;
% data2D = data2D_all(1:ndata,1:ndata) ;
% in_2D_data = data2D;

if fit_1d
     shifter = (IRF_1D_shift_range);
   normz  =zeros(length(IRF_1D_shift_range),1);
for df = 1:length(shifter)
    
    shift = IRF_1D_shift_range(df);    
    shifted_irf = IRF_shifter(nIRF,shift);
    
%% shrink
   
[shIRF,sh_time] = shrinker(shifted_irf,IRF_time,shrink_by);
[shData]        = data1D;

shIRF = shIRF/max(shIRF);  
shData(shData<0)=0;
            
%%  def basis and model Kernel

Tau1 = sh_time;
Tau2 = Tau1;

T1 = gen_basis(InitTime_T1,FinalTime_T1,L_basis,'log');
T2 = T1;

%%   make Kernel and SVD and convolute

K1 = gen_kernel(Tau1,T1,shIRF,2,L_basis);
K2 = K1;

[U1, S1, V1] = svds(K1, Singular_vals);
[U2, S2, V2] = svds(K2, Singular_vals);  
        %% 1D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
% using fake for loop just because we can collapse it for easy view

      input_data1D = shData;
      dm = max(max(input_data1D));
      input_data1D = input_data1D ./ dm;
      
     for nrep = 1:nrepeat
        
        iter_dat = input_data1D + randn(size(input_data1D))*noise_amplitude1d;
        
        [ILT_1D,Alpha1D,Fit1D] = ILT_core_st_V99(iter_dat,input_alpha,U1,S1,V1,1,1,1,cond_num_bound);
        
        data_all(:,nrep) = iter_dat;
        alpha_all(:,nrep) = Alpha1D.alpha;
        ILT_1D_all(:,nrep) = ILT_1D; %#ok<*SAGROW>
        Fit1D_all(:,nrep) = Fit1D;
        residues = iter_dat-Fit1D;
        res1_all(:,nrep) = residues;
        chi(nrep) = norm(residues)^2;
        normz(df) = mean(chi);
    end
        toc
        mean_1dILT = mean(ILT_1D_all,2);
        
        n_ilt = mean_1dILT'.*T1;
        n_ilt = n_ilt/sum(n_ilt);
%% plot 1D

    f4 = figure(6); 
    sgtitle(str_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
            semilogx(T1,n_ilt,'LineWidth',2)
            title('1DILT')
        subplot(222)
            semilogy(Tau1,input_data1D)
            title('input target')
            
        subplot(223)
            plot(res1_all)
            title('Residues')

        subplot(224)
            semilogy(sh_time,Fit1D_all,sh_time,abs(data_all),'.-',sh_time,shIRF)
            title('Fit and inputdata')
            
end   
drawnow
%% generate X2
    if shifter(1)~= shifter(end)
         figure(4)
        plot(shifter,normz,'.-')
        title('norm wrt IRF shifting')
        xlabel('IRF shift (ns)')
        ylabel('Chi sq')
        drawnow
    end
end
end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fit_2d
for MAIN_LOOP = 1:n_files
    
wt = waitbar(MAIN_LOOP/n_files);

cd(main_folder)

data2D = new_part{MAIN_LOOP}.sh2D;
% data2D_all = main_exp{MAIN_LOOP}.partition.TwoD_d1d1;
% data2D = data2D_all(1:ndata,1:ndata) ;

in_2D_data = process_data2D(data2D);

%% 2D analysis
 
%% shift IRF    
shift = IRF_2D_shift;
sh_irf = IRF_shifter(nIRF,shift);

%% shrink
[shIRF,sh_time] = shrinker(sh_irf,IRF_time,shrink_by);
shIRF = shIRF/max(shIRF);
       
%% data
        
shdata2D = in_2D_data;
shdata2D(shdata2D<0) = 0;

%% INIT
% make data
    Tau1 = sh_time;
    Tau2 = Tau1;

    T1 = gen_basis(InitTime_T1,FinalTime_T1,L_basis,'log');
    T2 = T1;      

    
%%   DEFINE KERNEL
  tic 
    K1 = gen_kernel(Tau1,T1,shIRF,2,L_basis);
    K2 = K1;

    [U1, S1, V1] = svds(K1, Singular_vals);
    [U2, S2, V2] = svds(K2, Singular_vals); 
        
        input_data2D = shdata2D;
        dm = max(max(input_data2D));
      input_data2D = input_data2D ./ dm;
 % %         

%% Robustness with noise       
%         
iter_dat = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
nilt = zeros(length(T1),length(T2),nrepeat);
res2_all = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
Fit2D_all = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
iter_alpha = zeros(nrepeat,1);

for kkk = 1:nrepeat
    
    kkk
    noice = randn(size(input_data2D));
    ampl = noise_amplitude2d*max(max(max(input_data2D)));
    iter_dat(:,:,kkk) = (ampl*noice) + input_data2D;
    iter_dat(iter_dat<0)=0;
    
        [FEst,Alpha_heel,Fit] = ILT_core_st_V99(iter_dat(:,:,kkk),input_alpha,U1,S1,V1,U2,S2,V2,cond_num_bound);
        
        iter_alpha(kkk) = Alpha_heel.alpha;
        ILT_2D_all(:,:,kkk) = FEst;
        Fit2D_all(:,:,kkk) = Fit;
        temp_res = iter_dat(:,:,kkk)-Fit;
        res2_all(:,:,kkk) = temp_res;
end
    
   mean_2dILT = mean(ILT_2D_all,3);
    log_mean_2dILT = log(mean_2dILT +1E-7);
    
      n2d = mean(ILT_2D_all,3).*T1.*T2';
    n2d = n2d/sum(sum(n2d));
    
    f3 = figure(5); 
    sgtitle(str_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
        cool_plot(T1,T2,log10(n2d + 1E-7))
        title ([str_name 'average fit'])
       
        subplot(222)
        mesh(iter_dat(:,:,kkk))
%         set(gca,'ZScale','log')
        title([str_name 'input data'])
%         
        subplot(223)
        mesh(res2_all(:,:,kkk))
%         set(gca,'ZScale','log')
        title([ str_name 'Residues'])
%             imagesc(fmodel)
        
        subplot(224)
        mesh(Fit2D_all(:,:,kkk))
%         set(gca,'ZScale','log')
        title([str_name 'Fit'])
%         mesh(rob_ilt)

%% save analyzed data
cd(dat_fold)
cd(new_dir)

% r = strfind(fname,'.');
% fname(r) = 'd';

parameters.T1 = T1;
parameters.T2 = T2;
parameters.Tau1 = Tau1;
parameters.Tau2 = Tau2;
parameters.microtime_resolution = microtime_resolution;
parameters.InitTime_T1 = InitTime_T1;
parameters.FinalTime_T1 = FinalTime_T1;                         
parameters.L_basis = L_basis;    
parameters.Singular_vals = Singular_vals;   
parameters.input_alpha = input_alpha; 
parameters.conv_alpha = mean(iter_alpha);
parameters.ndata = ndata;  
parameters.IRF_flag = IRF_flag;
parameters.nrepeat = nrepeat; 
parameters.noise_amplitude1d = noise_amplitude1d;
parameters.noise_amplitude2d = noise_amplitude2d;

batch_save.parameters = parameters;
if fit_1d
    batch_save.ilt1d = ILT_1D_all;
    batch_save.input1d = data1D;
    batch_save.res1d = res1_all;
    batch_save.fit1d = Fit1D_all;
end
if fit_2d
    batch_save.logilt2d = log_mean_2dILT;
    batch_save.meanilt2d = mean_2dILT;
    batch_save.input2d = data2D;
    batch_save.res2d = res2_all;
    batch_save.fit2d = Fit2D_all;
    ILT_frame(:,:,MAIN_LOOP) = n2d;
end
end
batch_save.fname = str_name;

sav_name = ['analyzed_' str_name];
save(sav_name ,'batch_save' )

% if fit_2d
% ImageJ;
frame_name = ['ILT_frames_' str_name];


% IJM.show('ILT_frame')
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
p26 = [100E-5 200E-3];
p27 = [200E-3 300E-3];
p28 = [300E-3 500E-3];
p29 = [500E-3 700E-3];
p30 = [700E-3 1];
p31 = [1 2];
p32 = [2 3];


frames = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20;p21;p22;p23;p24;p25;p26;p27;p28;p29;p30;p31;p32];
    dframe = frames(:,2) - frames(:,1);
dt = mean(frames,2);
for i = 1:32
pd(:,:,i) = ILT_frame(:,:,i)./dframe(i);
end

end
close(wt)

frame_name_pd = ['frame_pd_' str_name];
save(frame_name,'ILT_frame','dt','frames')


% IJM.show('ILT_frame')
imshow3D(ILT_frame)
fprintf('***************** ILT batch analysis complete ***************** \n')

%% functions 

function [T] = gen_basis(init_time,final_time,L_basis,basis_scale)
    switch basis_scale
           case 'log'
                T = logspace(log10(init_time), log10(final_time), L_basis);
           case 'linear'
                 T = linspace((init_time),(final_time), L_basis);
           otherwise
           error('ERROR: specify time scale properly')
     end
end

function [K] = gen_kernel(Tau,T,in_IRF,IRF_flag,L_basis)

    global microtime_resolution
    ndata = length(Tau);
 Kernel = @(time,TimeConst)(((1./TimeConst).*(exp(- time * (1./ TimeConst)))));
 
    switch IRF_flag
        
        case 0 
            L = Kernel (Tau,T);         
            
        case 1
            resolution = microtime_resolution;
            bckg = 0;
            power2 = 7E0;
            IRF_width2 = 2*resolution;
            IRF_Amp2 = 1.665*IRF_width2;
            IRF_t02 = 30*resolution;
            IRF_t2 = (1:256)*resolution; 
            tau2 = 0.05;

            term1_const2 = (sqrt(pi)*IRF_width2*IRF_Amp2/2);
            term2_ex2  = exp(((IRF_t02-IRF_t2)/tau2) + (IRF_width2.^2)/4*tau2.^2);
            term3_eror2 = (erfc(((IRF_t02-IRF_t2)/IRF_width2) +IRF_width2/2*tau2));

            ck2 = power2*term1_const2.*(term2_ex2) .* term3_eror2 +bckg;

            sim_irf = ck2'/max(ck2);
            % plot(sim_irf)  

            C1 = Kernel (Tau,T);
            TK1 = zeros(ndata,L_basis);

            for con = 1:L_basis
                temp_conv1 = conv(C1(:,con),sim_irf);
                TK1(:,con) = temp_conv1(1:ndata);

            end
                K = TK1;
                
        case 2 % custom experimental input IRF
            C1 = Kernel (Tau,T);
            TK1 = zeros(length(Tau),L_basis);
   
            for con = 1:L_basis
                temp_conv1 = conv(C1(:,con),in_IRF);
                TK1(:,con) = temp_conv1(1:ndata);
            end
            K = TK1;
                
    end
end

function [nIRF,IRF_time] = load_IRF(IRF_file)
    global microtime_resolution
    global ndata
    
%         IRF_file = 'irf_wavelength_p3_40mhz.mat';
        load(IRF_file)
        IRF_selected = 'd1w500';
        in_irf = d1w500(1:4093);
        in_irf = in_irf - mean(in_irf(1:250));
        
        IRF_time = (1:ndata)*microtime_resolution;
        nIRF = in_irf/max(in_irf);
%         nIRF(nIRF<0) = 0;
end

function [sh_IRF] = IRF_shifter(in_IRF,shift)
    % shifts IRF horizontally
    
    if sign(shift)>=0
        sh_IRF = [zeros(shift,1); in_IRF(1:end-shift)];
    elseif sign(shift) == -1
        sh_IRF = [(in_IRF(abs(shift-1):end)) ; zeros(abs(shift),1)];
   end
end

function [sh_trace,sh_time] = shrinker(in_trace,in_time,shrink_by)
    % data binning by a specified factor
    
    dat_size = length(in_trace);
    shrinked_size = ceil(dat_size/shrink_by);
    shrink = 1:shrink_by:dat_size;
    
    sh_trace = zeros(shrinked_size,1);
    sh_time = zeros(shrinked_size,1);
    for i =1:shrinked_size-1
        sh_trace(i) = sum(in_trace(shrink(i):(shrink(i+1)-1)));
        sh_time(i) = in_time(shrink(i));
    end
    sh_trace(shrinked_size) = sum(in_trace(shrink(shrinked_size):length(in_trace)));
    sh_time(shrinked_size) = in_time(end);    
end

function [proc_data] = process_data(data1D)

    data1D = data1D-mean(data1D(1:10));
    data1D(160:end) = 0;
     data1D = data1D';

    edit_bin_edges = 1;
    if edit_bin_edges
        data1D(1) =[];
        data1D(end) = data1D(end)/2;
        data1D(end+1) = data1D(end);
    end

    if data1D(1)==0
        error('ERROR: Possible miss calculated data_bins; turn ON edit_bin_edges');
    end
    proc_data = data1D;
end

function [pdata2D] = process_data2D(data2D)
    
    data2D(160:256,:) = 0;
    data2D(:,160:256) = 0;
    pdata2D = data2D;
    
end