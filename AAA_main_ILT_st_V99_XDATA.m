% figure out monte analysis first
clc
clear                                                       % Clear all variables
tic
%% USER DEFINED

t_name = 'EGFP_glycerol'
   fit_1d = 1;
   fit_2d = 0;

IRF_1D_shift_range = -41;
IRF_2D_shift = -35;

%% init

microtime_resolution = 0.004;
% base_microtime_resolution = 0.004;
   shrink_by = 1;

InitTime_T1 = 0.01;
     FinalTime_T1 = 50;                         
     L_basis = 200;
     
     Singular_vals = 40;
     cond_num_bound = 1E9;
     input_alpha = 0;

   ndata = 4093;

%    IRF_flag = 1;
   fixed_basis_flag = 0;
   if  fixed_basis_flag
       L_basis = 3;
   end
   
   nrepeat = 1; 
   noise_amplitude1d =0E-3;
   noise_amplitude2d = 0E-4;
   
       basis_scale = 'log';
       
%% code
for code_dummy =1
current_folder = pwd; 

cd('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\Lee_data\Lee_EGFP_glycerol')

[select_file, file_path] = uigetfile;
load([file_path select_file])
% load('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\20201019_rfp_lee_mini\20201019_rfp_lee_mini-42_1_2_3.mat ')
cd(current_folder);

%% IRF

IRF_file = 'irf_wavelength_p3_40mhz.mat';
    load(IRF_file)
    IRF_selected = 'd2w488';
    in_irf1 = d1w488(1:4093);
    in_irf2 = in_irf1;

    
% data1D = partition.OneD_dec1 + partition.OneD_dec2 ;
 data1D = partition.OneD_dec1;
edit_bin_edges = 1;
if edit_bin_edges
    data1D(1) =[];
    data1D(end) = data1D(end)/2;
    data1D(end+1) = data1D(end);
end

if data1D(1)==0
    error('ERROR: Possible miss calculated data_bins; turn ON edit_bin_edges');
end

% data2D = partition.TwoD_d1d1 + partition.TwoD_d1d2+partition.TwoD_d2d1+partition.TwoD_d2d2;
data2D = partition.TwoD_d1d1;


   %% Prep
   
     IRF_time = (1:ndata)*microtime_resolution;       
%           ex_IRF1 = in_irf1-mean(in_irf1(100:300));
          ex_IRF1 = in_irf1;
          nirf1 = ex_IRF1/max(ex_IRF1);
          nirf1(nirf1<0) = 0;
          
%           ex_IRF2 = in_irf2-mean(in_irf2(100:300));
          ex_IRF2 = in_irf2;
          nirf2 = ex_IRF2/max(ex_IRF2);
          nirf2(nirf2<0) = 0;
   
          data1D = data1D - mean(data1D(100:300));
          in_2D_data = data2D;
end      
   %% 1D analysis using shifter loop 
   if fit_1d
   shifter = (IRF_1D_shift_range);
   normz  =zeros(length(shifter),1);
   shift_ilts = zeros(length(IRF_1D_shift_range),L_basis);
   for df = 1:length(shifter)
       
       shift1 = shifter(df);
       
   if sign(shift1)>=0
        IRF1 = [zeros(shift1,1); nirf1(1:end-shift1)];
    elseif sign(shift1) == -1
        IRF1 = [(nirf1(abs(shift1-1):end)) ; zeros(abs(shift1),1)];
   end
    
   %% shrink
   
           dat_size = length(nirf1);
            shrinked_size = ceil(dat_size/shrink_by);

            shrink = 1:shrink_by:dat_size;
            shirf1 = zeros(shrinked_size,1);
            shirt = zeros(shrinked_size,1);
            shpp2 = zeros(shrinked_size,1);
            
        for i =  1:shrinked_size-1
            shirf1(i) = sum(IRF1(shrink(i):(shrink(i+1)-1)));
            shirt(i) = IRF_time(shrink(i));
            shpp2(i) = sum(data1D(shrink(i):(shrink(i+1)-1)));
        end

        shirf1(shrinked_size) = sum(IRF1(shrink(shrinked_size):length(nirf1)));
        shirt(shrinked_size) = IRF_time(end);
        shpp2(shrinked_size) = sum(data1D(shrink(shrinked_size):length(nirf1)));

    shirf1 = shirf1/max(shirf1);  

    shdata= shpp2;
    shdata(shdata<0)=0;

            
            
    %%  def basis and model Kernel
%    shirt = (1:ndata(1))*microtime_resolution;
       ndata = shrinked_size; 

       Tau1 = shirt;
       Tau2 = Tau1;

   switch basis_scale
       case 'log'
            T1 = logspace(log10(InitTime_T1), log10(FinalTime_T1), L_basis);
       case 'linear'
             T1 = linspace((InitTime_T1),(FinalTime_T1), L_basis);
       otherwise
           error('ERROR: specify time scale properly')
   end
   if fixed_basis_flag
    T1 = [ 0.043 0.7565 2.477 ]
   end
   
   T2 = T1;
       Kernel_1 = @(Tau,TimeConst)((exp(- Tau * (1./ TimeConst))));

       %%   make Kernel and SVD and convolute

        C1 = Kernel_1 (Tau1,T1);
    C2 = Kernel_1(Tau2,T2);

        TK1 = zeros(ndata,L_basis);
        TK2 = zeros(ndata,L_basis);
   
        for con = 1:L_basis
            temp_conv1 = conv(C1(:,con),shirf1);
            TK1(:,con) = temp_conv1(1:ndata);

            temp_conv2 = conv(C1(:,con),shirf1);                       
            TK2(:,con) = temp_conv2(1:ndata);
        end
    
            K1 = TK1;
            K2 = TK2;

        [U1, S1, V1] = svds(K1, Singular_vals);
        [U2, S2, V2] = svds(K2, Singular_vals);    
 %% 1D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
% using fake for loop just because we can collapse it for easy view

if fit_1d
for data_1D_dummy = 1
      input_data1D = shdata';
      dm = max(max(input_data1D));
      input_data1D = input_data1D ./ dm;
      
     for nrep = 1:nrepeat
        
        iter_dat = input_data1D + randn(size(input_data1D))*noise_amplitude1d;
        
        [ILT_1D,Alpha1D,Fit1D,NoiseStd_1D] = ILT_core_st_V99(iter_dat,input_alpha,U1,S1,V1,1,1,1,cond_num_bound);
        
        data_all(:,nrep) = iter_dat;
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
        shift_ilts(df,:) = n_ilt;
%% plot 1D

    f4 = figure(17); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
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
            semilogy(shirt,Fit1D_all,shirt,abs(data_all),'.-',shirt,shirf1)
            ylim([1E-5 1.1])
            title('Fit and inputdata')
            
end   
drawnow
end
%% generate X2
    if length(shifter)>1
         figure(4)
        plot(shifter,normz,'.-')
        title('norm wrt IRF shifting')
        xlabel('IRF shift (ns)')
        ylabel('Chi sq')
        
        figure(10)
        plot(mean(shift_ilts),'.-')
        title('mean IRF shifted ILT')
    end
   end        
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D analysis
if fit_2d
for data_2D_dummy = 1
    
   shifter2 = (IRF_2D_shift);
   normz2D  =zeros(length(shifter2),1);
   shift_ilts2D = zeros(L_basis,L_basis,length(IRF_2D_shift));
for df2 = 1:length(shifter2)
%% shift   
shift1 = shifter2(df2);
shift2 = shifter2(df2);

if sign(shift1)>=0
    IRF1 = [zeros(shift1,1); nirf1(1:end-shift1)];
elseif sign(shift1) == -1
    IRF1 = [(nirf1(abs(shift1-1):end)) ; zeros(abs(shift1),1)];
end

if sign(shift2)>=0
    IRF2 = [zeros(shift2,1); nirf2(1:end-shift2)];
elseif sign(shift2) == -1
    IRF2 = [(nirf2(abs(shift2-1):end)) ; zeros(abs(shift2),1)];
end
%% shrink
    dat_size = length(nirf1);
        shrinked_size = ceil(dat_size/shrink_by);

        shrink = 1:shrink_by:dat_size;
        shirf1 = zeros(shrinked_size,1);
        shirf2 = zeros(shrinked_size,1);
        shirt = zeros(shrinked_size,1);
%         shpp2 = zeros(shrinked_size,1);   
        

            for i =  1:shrinked_size-1
                shirf1(i) = sum(IRF1(shrink(i):(shrink(i+1)-1)));
                shirf2(i) = sum(IRF2(shrink(i):(shrink(i+1)-1)));
                shirt(i) = IRF_time(shrink(i));
%                 shpp2(i) = sum(dec_sum(shrink(i):(shrink(i+1)-1)));
            end

            shirf1(shrinked_size) = sum(IRF1(shrink(shrinked_size):length(IRF1)));
            shirf2(shrinked_size) = sum(IRF2(shrink(shrinked_size):length(IRF2)));
            shirt(shrinked_size) = IRF_time(end);
%             shpp2(shrinked_size) = sum(dec_sum(shrink(shrinked_size):length(nirf1)));

        shirf1 = shirf1/max(shirf1);
        shirf2 = shirf2/max(shirf2);
        
%          shpp2 = shpp2/max(shpp2);
%             shdata = shpp2*shpp2';
       
        %% data
        
        dat_size = length(nirf1);
        dat_shrink = (1:shrink_by:dat_size);
        shrinked_size = ceil(dat_size/shrink_by);  
        
    shrink_datvec = zeros(shrinked_size,shrinked_size);
%     unc_dat = zeros(shrinked_size,shrinked_size);
    
    for i = 1:shrinked_size-1
        for j = 1:shrinked_size-1
            shrink_datvec(i,j) = sum(sum(in_2D_data(shrink(i):(shrink(i+1)-1),shrink(j):(shrink(j+1)-1))));
%             unc_dat(i,j) = sum(sum(unc_b1(shrink(i):(shrink(i+1)-1),shrink(j):(shrink(j+1)-1))));
        end
    end      
    
         shdata = shrink_datvec;
%          shdata = shrink_datvec/max((max(shrink_datvec)));  
        shdata(shdata<0) = 0;

%% INIT
% make data
          
    ndata = [ shrinked_size shrinked_size]; % number of data points
          Tau1 = shirt;
          Tau2 = Tau1;
    
   switch basis_scale
       case 'log'
            T1 = logspace(log10(InitTime_T1), log10(FinalTime_T1), L_basis);
       case 'linear'
             T1 = linspace((InitTime_T1),(FinalTime_T1), L_basis);
       otherwise
           error('ERROR: specify time scale properly')
   end
     if fixed_basis_flag
       T1 = [ 0.043 0.7565 2.477 ]
     end
   T2 = T1;
   
   Kernel_1 = inline('exp(- Tau * (1./ TimeConst))','Tau','TimeConst');
%%   DEFINE KERNEL
  tic 
  
	C_1 = Kernel_1 (Tau1,T1);
    
    TK1 = zeros(shrinked_size,L_basis);
    TK2 = zeros(shrinked_size,L_basis);
    for con = 1:L_basis
        temp_conv1 = conv(C_1(:,con),shirf1);
        TK1(:,con) = temp_conv1(1:ndata(1));
        
        temp_conv2 = conv(C_1(:,con),shirf2);
        TK2(:,con) = temp_conv2(1:ndata(2));
    end

    %IRF convoluted kernel
    K1 = TK1;
    K2 = TK2;

        [U1, S1, V1] = svds(K1, Singular_vals);
        [U2, S2, V2] = svds(K2, Singular_vals);
        
 % %         

%% Robustness with noise       
%         
        input_data2D = shdata;
      dm = max(max(input_data2D));
      input_data2D = input_data2D ./ dm;

iter_dat = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
nilt = zeros(length(T1),length(T2),nrepeat);
iter_res = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
iter_fit = zeros(size(input_data2D,1),size(input_data2D,2),nrepeat);
iter_alpha = zeros(nrepeat,1);

for kkk = 1:nrepeat
    
    kkk
   iter_dat(:,:,kkk) = input_data2D + randn(size(input_data2D))* noise_amplitude2d;
%     iter_dat(iter_dat<0)=0;
    
        [FEst,Alpha_heel,Fit,NoiseStd_2D] = ILT_core_st_V99(iter_dat(:,:,kkk),input_alpha,U1,S1,V1,U2,S2,V2,cond_num_bound);
        
        iter_alpha = Alpha_heel.alpha;
        ILT_2D_all(:,:,kkk) = FEst;
        iter_fit(:,:,kkk) = Fit;
        temp_res = iter_dat(:,:,kkk)-Fit;
        iter_res(:,:,kkk) = temp_res;
        chi(:,kkk) = norm(temp_res).^2;
end
    normz2D(df2) = mean(chi);

   mean_2dILT = mean(ILT_2D_all,3);
    log_mean_2dILT = log(mean_2dILT +1E-7);
    
      n2d = mean(ILT_2D_all,3).*T1.*T2';
    n2d = n2d/sum(sum(n2d));
    
    shift_ilts2D(:,:,df2) = n2d;
    
    f3 = figure(5); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
        cool_plot(T1,T2,n2d)
        title ([t_name 'average fit'])
       
        subplot(222)
        mesh(iter_dat(:,:,kkk))
%         set(gca,'ZScale','log')
        title([t_name 'input data'])
%         
        subplot(223)
        mesh(iter_res(:,:,kkk))
%         set(gca,'ZScale','log')
        title([ t_name 'Residues'])
%             imagesc(fmodel)
        
        subplot(224)
        mesh(iter_fit(:,:,kkk))
%         set(gca,'ZScale','log')
        title([t_name 'Fit'])
%         mesh(rob_ilt)
end 
end
end