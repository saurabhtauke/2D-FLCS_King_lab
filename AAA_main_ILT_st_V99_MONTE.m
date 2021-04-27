clc
clearvars -except IJM                                                      % Clear all variables
tic

t_name = 'monte-'
% t_name = 'Monte-IRF-t0-shift-1-4'

if isfile(['analyzed_' t_name '.mat'])
%     error( ' ERROR : file exists')
end

  ndata = 250;
  
current_folder = pwd; 
cd('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\')

[select_file, file_path] = uigetfile;
load([file_path select_file])
% load('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\Resolution_2_peaks\resolution-1-3\resolution-1-322_1_0.1_0.5.mat')
% load('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\Resolution_IRF_2Peaks\IRF-25s_1-3\IRF-25s_1-3-22_1_0.1_0.5.mat')
cd(current_folder)

data1D = partition.OneD_dec1 + partition.OneD_dec2 ;
edit_bin_edges = 1;
if edit_bin_edges
    data1D(1) =[];
    data1D(end) = data1D(end)/2;
    data1D(end+1) = data1D(end);
end

data1D = data1D(1:ndata);
% if data1D(1)==0
%     error('ERROR: Possible miss calculated data_bins; turn ON edit_bin_edges');
% end

data2D_all = partition.TwoD_d1d1 + partition.TwoD_d1d2+partition.TwoD_d2d1+partition.TwoD_d2d2;
data2D = data2D_all(1:ndata,1:ndata) ;
%% init

microtime_resolution = 0.064;

InitTime_T1 = 0.05;
     FinalTime_T1 =15;                         
     L_basis = 100;
     
     Singular_vals = 30;
      cond_num_bound = 1E16; % 1E6 is ok
      
     input_alpha = 0    
   
   fit_1d = 1;
   fit_2d = 1;
  
   nrepeat = 1; 
   noise_amplitude1d = 0E-3;
   noise_amplitude2d = 0E-3;
  
   IRF_flag = 0;
   itrable = 30+(-2:2)*2*0.15;
   IRF_1D_shift_range = itrable ;
   IRF_t0 = 0.05;
   
   loop_range  = 1%length(itrable);
   
   basis_scale = 'log';
   
   all_ilt = zeros(loop_range,L_basis);
   %%  def basis and model Kernel
   shirt = (1:ndata(1))'*microtime_resolution;
   
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
%    T1 = [0.2 0.8 2 4.5];
   T2 = T1;
   
   Kernel_1 = @(Tau,TimeConst)((exp(- Tau * (1./ TimeConst))));
   
   for multi_iter = 1:loop_range
       IRF_1D_shift_range = itrable(multi_iter);
       % change any variable we want to loop over in this section
       
%% IRF synth
switch IRF_flag
    case 1
        resolution = microtime_resolution;
        bckg = 0;
        power2 = 7E0;
        IRF_width2 = 2*resolution;
        IRF_Amp2 = 1.665*IRF_width2;
        IRF_t02 = IRF_1D_shift_range*resolution;
        IRF_t2 = (1:256)*resolution; 
        tau2 = IRF_t0;

        term1_const2 = (sqrt(pi)*IRF_width2*IRF_Amp2/2);
        term2_ex2  = exp(((IRF_t02-IRF_t2)/tau2) + (IRF_width2.^2)/4*tau2.^2);
        term3_eror2 = (erfc(((IRF_t02-IRF_t2)/IRF_width2) +IRF_width2/2*tau2));

        ck2 = power2*term1_const2.*(term2_ex2) .* term3_eror2 +bckg;

        sim_irf = ck2'/max(ck2);
        % plot(sim_irf)  

%%      make Kernel and SVD and convolute
        
        C1 = Kernel_1 (Tau1,T1);
        C2 = Kernel_1(Tau2,T2);

        TK1 = zeros(ndata,L_basis);
        TK2 = zeros(ndata,L_basis);
   
        for con = 1:L_basis
            temp_conv1 = conv(C1(:,con),sim_irf);
            TK1(:,con) = temp_conv1(1:ndata);

            temp_conv2 = conv(C1(:,con),sim_irf);                       
            TK2(:,con) = temp_conv2(1:ndata);
        end
    
            K1 = TK1;
            K2 = TK2;

    case 0 
        C1 = Kernel_1 (Tau1,T1);
        C2 = Kernel_1(Tau2,T2);
        
            K1 = C1;
            K2 = C2;
end

        [U1, S1, V1] = svds(K1, Singular_vals);
        [U2, S2, V2] = svds(K2, Singular_vals);
        
%% 1D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
% using fake for loop just because we can collapse it for easy view

if fit_1d
   shifter = (IRF_1D_shift_range);
   normz  =zeros(length(shifter),1);
   
for model_1D_dummy = 1
   
%    for df = 1:length(shifter)
%        
%        shift1 = shifter(df);
%        
%     if sign(shift1)>=0
%         IRF1 = [zeros(shift1,1); nirf1(1:end-shift1)];
%     elseif sign(shift1) == -1
%         IRF1 = [(nirf1(abs(shift1-1):end)) ; zeros(abs(shift1),1)];
%     end
    
      input_data1D = data1D;
      dm = max(max(input_data1D));
      input_data1D = input_data1D ./ dm;

      
%% fitting 1D

    for nrep = 1:nrepeat
        
        iter_dat = input_data1D + randn(size(input_data1D))*(noise_amplitude1d);
        
         mtilde = (U1' * iter_dat' * 1);
         mtilde2 = U1 *mtilde* 1';
         Noise_estimate1D(multi_iter) = mean(std(mtilde2 - iter_dat'))
        
        [ILT_1D,Alpha1D,Fit1D,NoiseStd_1D] = ILT_core_st_V99(iter_dat,input_alpha,U1,S1,V1,1,1,1,cond_num_bound);
        
        data_all(:,nrep) = iter_dat;
        ILT_1D_all(:,nrep) = ILT_1D; %#ok<*SAGROW>
        Fit1D_all(:,nrep) = Fit1D;
        res1_all(:,nrep) = iter_dat-Fit1D;
         alpha_all(:,nrep) = Alpha1D.alpha;
        
        chi(multi_iter) = norm(res1_all)^2;
    end
        toc
        mean_1dILT = mean(ILT_1D_all,2);
        
        n_ilt = mean_1dILT'.*T1;
%         n_ilt = n_ilt/sum(n_ilt);
        
    all_ilt(multi_iter,:) = n_ilt;

%% plot 1D

    f4 = figure(8); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
%             semilogx(T1,mean_1dILT,'-o',T1,ILT_1D_all)
            semilogx(T1,n_ilt,'LineWidth',2)
            title(['1DILT \alpha = '  num2str(mean(alpha_all),'%10.5e\n')])
%             legend([t_name '-shift-' num2str(IRF_1D_shift_range)])
%             xlim([0.8 7])
        subplot(222)
            semilogy(Tau1,iter_dat)
            title('input target')
            
        subplot(223)
            plot(res1_all)
%             ylim([-0.015 0.015])
            title('Residues')

        subplot(224)
            semilogy(shirt,abs(data_all),'.-',shirt,Fit1D_all)
%             ylim([1E-4 1])
            title('Fit and inputdata')
            
end   
drawnow
end

%% 2D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if fit_2d
for model_2D_dummy = 1

      
      input_data2D = data2D;
      dm = max(max(input_data2D));
      input_data2D = input_data2D ./ dm;
      
%% fitting 2D
    for nrep2 = 1:nrepeat
        data_tofit2D = input_data2D + randn(size(input_data2D))* noise_amplitude2d;

        [ILT_2D,Alpha,Fit,NoiseStd_2D] = ILT_core_st_V99(data_tofit2D,input_alpha,U1,S1,V1,U2,S2,V2,cond_num_bound);
        
        data2_all(:,:,nrep2) = data_tofit2D;
        ILT_2D_all(:,:,nrep2) = ILT_2D; %#ok<*SAGROW>
        Fit2D_all(:,:,nrep2) = Fit;
        res2_all(:,:,nrep2) = data_tofit2D-Fit;
        alpha_all2d(:,nrep2) = Alpha.alpha;
       
        chi(multi_iter) = norm(res2_all(:,:,nrep2))^2;
    end
    
    mean_2dILT = mean(ILT_2D_all,3)+1E-7;
    
    n2d = mean(ILT_2D_all,3).*T1.*T2';
    n2d = (n2d/sum(sum(n2d)));
        toc
%% plot 2D

    f5 = figure(15); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
            cool_plot(T1,T2,n2d)
%                 imagesc(T1,T2,mean_2dILT)
            title('log(2DILT)')
        subplot(222)
            mesh((Tau1),(Tau2),data_tofit2D)
            title('input target')
    %         
        subplot(223)
            mesh((Tau1),(Tau2),res2_all(:,:,1))
            title('Residues')

        subplot(224)
            mesh((Tau1),(Tau2),Fit2D_all(:,:,1))
            title('Fit')
end
end
%%
   end
   
   if multi_iter > 1 
       figure
       hold on
       loglog(itrable,chi,'.-','LineWidth',1.5)
   end