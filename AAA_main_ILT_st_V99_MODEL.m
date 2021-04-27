clc
clearvars -except IJM   % Clear all variables

t_name = 'MODEL'

current_folder = pwd; 
cd('D:\OneDrive - unist.ac.kr\work\ILT\methods_v_inf_final\preprocess\Partition\saurabh_data')
[select_file, file_path] = uigetfile;
load([file_path select_file])
cd(current_folder)


x = partition.OneD_dec1+partition.OneD_dec2;
x2 = partition.TwoD_d1d1+partition.TwoD_d2d1+partition.TwoD_d1d2+partition.TwoD_d2d2;

%% init

  for ttt = 1:1  
%    rng(30)   
microtime_resolution = 0.064;
 
InitTime_T1 = 0.1;
     FinalTime_T1 = 10;                         
     L_basis = 100;
      
     Singular_vals = 30;   
     cond_num_bound = 1E9; % 1E6 is ok
     
     input_alpha = 0;  
% input_alpha = itrable(ttt)

   ndata = 256;
   
   fit_1d = 0;
   fit_2d = 1;
   
   ana_1d = 0;
   ana_2d = 0;

    IRF_flag = 0;
   
   nrepeat = 1;
  
   noise_amplitude1d = 1E-2;
   noise_amplitude2d = 1E-2;
   
   basis_scale = 'log'

   %%  def basis and model
   data_time = (1:ndata)*microtime_resolution;
   
   Tau1 = data_time';
   Tau2 = Tau1;
   
   T1 = gen_basis(InitTime_T1,FinalTime_T1,L_basis,'log');
   
   T2 = T1;
   
%%   make Kernel and SVD
  
	
K = gen_kernel(Tau1,T1,IRF_flag,microtime_resolution);
% K1 = 2*(sqrt(K+(3/8)));
% K1 = sqrt(K);
K1 = K;
K2 = K1;

        [U1, S1, V1] = svds(K1, Singular_vals);
        [U2, S2, V2] = svds(K2, Singular_vals);
        
        
%% 1D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
% using fake for loop just because we can collapse it for easy view

if fit_1d
for model_1D_dummy = 1
%% make model 1D
    %~~~~~~~~~~1Dmodel~~~~~~~~~~~~~~~~~~~~~~~~~
    fmodel1D = ones(1,L_basis)*1E-5;
    amplitude = 1e3
    
%         fmodel1D(50) = 1;
%         fmodel1D(70) = 1;
%         fmodel1D(30) = 1;

fmodel1D = exp(-([-39:60]./2 ).^2/2) + exp(-([-69:30]./2 ).^2/2);
% fmodel1D = exp(-([-70:29]./2 ).^2/2) + exp(-([-20:79]./2 ).^2/2) + exp(-([-45:54]./2 ).^2/2);
          fmodel1D = amplitude*fmodel1D'./sum(fmodel1D);
%           
    data1D = (K)*fmodel1D;  
    
    input_data1D = poissrnd(data1D');
    dm = max(input_data1D);
   
input_data1D = x;

%% fitting 1D

    for nrep = 1:nrepeat
        
%         iter_dat = 2*sqrt((input_data1D)+(3/8));
%         iter_dat = sqrt(input_data1D);
        iter_dat = input_data1D;
        
        [ILT_1D,Alpha1D,Fit1D,NoiseStd_1D,oppar] = ILT_core_st_V99(iter_dat,input_alpha,U1,S1,V1,1,1,1,cond_num_bound);
        
        data_all(:,nrep) = iter_dat;
        ILT_1D_all(:,nrep) = ILT_1D; %#ok<*SAGROW>
        Fit1D_all(:,nrep) = Fit1D;
        res1_all(:,nrep) = iter_dat-Fit1D;
        norm_res(nrep) = norm(iter_dat-Fit1D);
        
    end
        
        mean_1dILT = mean(ILT_1D_all,2);
        n_ilt = mean_1dILT';

        
%% plot 1D

    f4 = figure(4); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
            plot(n_ilt,'.-')
            hold on
%             plot(fmodel1D,'k-o')
            title(['1DILT \alpha = '  num2str(oppar.apar(end),'%10.5e\n')])
        subplot(222)
            semilogx(T1,fmodel1D)
            title('input target')
    %         
        subplot(223)
            plot(res1_all)
%             ylim([-0.015 0.015])
            title('Residues')
            hold on

        subplot(224)
            semilogy(data_time,(data_all),'.-',data_time,Fit1D_all)
%             ylim([1E-4 1])
            title('Fit and inputdata')
            hold on
            
end
drawnow
end
%% 2D ILT analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if fit_2d
for model_2D_dummy = 1
%% make model 2D
    %~~~~~~~~~ 2Dmodel~~~~~~~~~~~~~~~~~~~~~~~~~~
    fmodel2D = ones(L_basis,L_basis)*1E-5; 
amp2 = 1E4;
    
%     fmodel1 = exp(-([-70:29]./2 ).^2/2) + exp(-([-20:79]./2 ).^2/2) + exp(-([-45:54]./2 ).^2/2);
%     fmodel2 = exp(-([-70:29]./2 ).^2/2) + exp(-([-20:79]./2 ).^2/2) + exp(-([-45:54]./2 ).^2/2);     
%     fmodel2D = fmodel1' * fmodel2;
    
%     fmodel2D = fmodel1D*fmodel1D';
%     fmodel2D = amp2*fmodel2D/sum(sum(fmodel2D));
%     
%     p1 = 30;
%     p2 = 57;
%     p3 = 75;
%     
%       fmodel2D(p1,p1) = amp2;
%         fmodel2D(p2,p2) = amp2;
%         fmodel2D(p3,p3) = amp2;
%         
%             fmodel2D(p1,p2) = amp2/5;
%             fmodel2D(p2,p3) = amp2/5;
%             fmodel2D(p3,p1) = amp2/5;
%             
%              fmodel2D(p2,p1) = amp2;
%             fmodel2D(p3,p2) = amp2;
%             fmodel2D(p1,p3) = amp2;

        
%       data2D = K1 * fmodel2D * K2'; 

        data2D = x2;

      input_data2D = data2D;
      
      
%% fitting 2D

    for nrep2 = 1:nrepeat
          
          data_tofit2D = poissrnd(input_data2D);
          iter2 = data_tofit2D;
        
        [ILT_2D,Alpha,Fit,NoiseStd_2D,oppar] = ILT_core_st_V99(data_tofit2D,input_alpha,U1,S1,V1,U2,S2,V2,cond_num_bound);
        
        data2_all(:,:,nrep2) = data_tofit2D;
        ILT_2D_all(:,:,nrep2) = ILT_2D; %#ok<*SAGROW>
        Fit2D_all(:,:,nrep2) = Fit;
        res2_all(:,:,nrep2) = data_tofit2D-Fit;
        alpha_all2d(:,nrep2) = Alpha.alpha;
       
    end
    toc
    mean_2dILT = mean(ILT_2D_all,3);
     n2d = mean(ILT_2D_all,3);
    n2d = n2d/sum(sum(n2d));
    

%% plot 2D

    f5 = figure(5); 
    sgtitle(t_name,'FontSize',29,'FontWeight','bold')
        subplot(221)
%             cool_plot(T1,T2,n2d)
              imagesc(mean_2dILT)
            title('mean_2DILT')
        subplot(222)
            mesh(Tau1,Tau2,data_tofit2D)
            title('input target')
    %         
        subplot(223)
            mesh(res2_all(:,:,1))
            title('Residues')

        subplot(224)
            mesh(Fit2D_all(:,:,1))
            title('Fit')
end
end
%%
stopper = oppar.stopper
c_len = size(oppar.cv_dump,2);

if ana_1d == 1
    
        evol_spec(:,:,:) = (max(input_data1D)).* (max(0,oppar.kk'*(oppar.cv_dump(:,:))));
    
    f26 = figure(26);
    subplot(3,4,1)
    mesh(mean(evol_spec,3))
    view([120 30])

    cchi = sum(abs((iter_dat'-K1*mean(evol_spec,3)))) ;

%     f27 = figure(27)
    subplot(3,4,2)
    loglog(oppar.apar,(cchi),'.-')
    title('log10(L1-distance)')
    hold on
    
    chi = @(dat,K,spec)(((dat-K*spec).^2));

    log_L = @(dat,K,spec)(0.5.*log(dat)+0.5*(chi(dat,K,spec)./(dat)));

    test_l = log_L(iter_dat',K1,evol_spec(:,1));

    loglik = zeros(length(iter_dat),size(evol_spec,2));

    for i = 1:size(evol_spec,2)
        loglik(:,i) = log_L(iter_dat',K1,evol_spec(:,i));
        loglik(isinf(loglik) | isnan(loglik))= 0;
        
        chisq(i) = sum(chi(iter_dat',K1,evol_spec(:,i)));
    end
    
    subplot(3,4,8)
    loglog(oppar.apar,chisq,'.-');title('\chi ^2');hold on
    
    sll = sum(loglik);
    % f28 = figure(28)
    subplot(3,4,3)
    loglog(oppar.apar,sll,'o-');title('neg-log - likelihood tikhonov');hold on

    et= -evol_spec.*log(evol_spec);
    et(isnan(et))=0;
    sumet = sum(et);

    % f29 = figure(29)
    subplot(3,4,4)
    semilogx(oppar.apar,sumet/sum(iter_dat),'.-');title('spec entropy');hold on
    
    subplot(3,4,11)
    loglog(oppar.apar,oppar.apar,'k.-');title('int-chi');hold on
    loglog(oppar.apar,oppar.fronorm.*oppar.apar/100,'.-')

    subplot(3,4,12)
    loglog(oppar.apar,oppar.fronorm,'.-');title('norm(c-vec)');hold on

    subplot(3,4,7)
    loglog(oppar.apar,oppar.fval,'.-');title('fval');hold on
    % subplot(3,4,8)

    % kl = @(f,g)(sum(f.*log(f./g)));
    % irad=@(f,g)(kl(f,(f+g)/2) + kl(g,(f+g)/2));

    kls1 = kl(evol_spec,fmodel1D);
    kls2 = kl(fmodel1D,evol_spec);
    kld1 = kl(iter_dat',K1*evol_spec);
    kld2 = kl(K1*evol_spec,iter_dat');

    irads = irad(evol_spec,fmodel1D);
    iradd = irad(iter_dat',K1*evol_spec);

    % figure(51)
    subplot(3,4,5)
    loglog(oppar.apar,kld1);title('kld1');hold on
    subplot(3,4,6)
    loglog(oppar.apar,kld2);title('kld2');hold on

    % figure(52)
    subplot(3,4,9)
    loglog(oppar.apar,irads,'.-');title('irads');hold on
    subplot(3,4,10)
    loglog(oppar.apar,iradd,'.-');title('iradd');hold on
    % loglog(oppar.apar,oppar.apar,'.')


    mt= oppar.mtilde2.*iter_dat(1);
    mt(mt<0)= 0 ;
    j1 = irad(iter_dat',mt);

    

end

if ana_2d ==1
    evol_spec2 = max(max(data2D)).*(max(0,reshape(oppar.kk'*(oppar.cv_dump(:,:)),L_basis,L_basis,c_len)));
    
    fit2 = zeros(length(data_time),length(data_time),c_len);
    cchi2 = zeros(1,c_len);
    
    for i  = 1:c_len
        fit2(:,:,i) = K2*evol_spec2(:,:,i)*K1';
        cchi2(i) = sum(sum((fit2(:,:,i) - data_tofit2D).^2));
    end
       
    figure(31)
    subplot(3,4,2)
    loglog(oppar.apar,cchi2,'.-');title('log10(L1-distance)');hold on

    chi2 = @(dat,K,spec)((dat-K*spec*K').^2);
    logl2 = @(dat,K,spec)(0.5.*log(dat)+0.5*(chi2(dat,K,spec)./(dat)));
    loglik2 = zeros([size(iter2),size(evol_spec2,3)]);
    
    for i = 1:size(evol_spec2,3)
        loglik2(:,:,i) = logl2(iter2,K1,evol_spec2(:,:,i));
        loglik2(isinf(loglik2) | isnan(loglik2))= 0;
        
         chisq2(i) = sum(sum(chi2(iter2,K1,evol_spec2(:,:,i))));
    end
    
    subplot(3,4,8)
    loglog(oppar.apar,chisq2,'.-');title('\chi ^2');hold on
    
    sll2 = sum(sum(loglik2));
    sll2 = squeeze(sll2);
    subplot(3,4,3)
    loglog(oppar.apar,sll2,'o-');title('neg-log - likelihood tikhonov');hold on
    
    et= -evol_spec2.*log(evol_spec2);
    et(isnan(et))=0;
    sumet = squeeze(sum(sum(et)));
        
    subplot(3,4,4)
    semilogx(oppar.apar,sumet/sum(sum(iter2)),'.-');title('spec entropy');hold on
    
    
    subplot(3,4,11)
    loglog(oppar.apar,oppar.apar,'k.-');title('int-chi');hold on
    loglog(oppar.apar,oppar.fronorm.*oppar.apar/100,'.-')

    subplot(3,4,12)
    loglog(oppar.apar,oppar.fronorm,'.-');title('norm(c-vec)');hold on

    subplot(3,4,7)
    loglog(oppar.apar,oppar.fval,'.-');title('fval');hold on
    
    kls1 = squeeze(sum(kl(evol_spec2,fmodel2D)));
    kls2 = squeeze(sum(kl(fmodel2D,evol_spec2)));
    kld1 = squeeze(sum(kl(iter2,fit2)));
    kld2 = squeeze(sum(kl(fit2,iter2)));

    irads = squeeze(sum(irad(evol_spec2,fmodel2D)));
    iradd = squeeze(sum(irad(iter2,fit2)));
    
    subplot(3,4,5)
    loglog(oppar.apar,kld1);title('kld1');hold on
    subplot(3,4,6)
    loglog(oppar.apar,kld2);title('kld2');hold on
    
    subplot(3,4,9)
    loglog(oppar.apar,irads,'.-');title('irads');hold on
    subplot(3,4,10)
    loglog(oppar.apar,iradd,'.-');title('iradd');hold on
    
    mt= oppar.mtilde2.*iter2(1);
    mt(mt<0)= 0 ;
    j1 = irad(iter2,mt);
 

    
%     evol_spec2 = 0;
    drawnow
end

% ns(ttt) = oppar.NoiseStd;
% counts(ttt) = sum(input_data1D);
% apar(ttt) = Alpha1D.alpha;
% ress(ttt) = sum(res1_all.^2)/256;

% lsq(ttt) = sum(sum((U1*S1*V1' - K1).^2));
% time(ttt) = oppar.time

% mt(:,ttt) = oppar.mtilde2;
% it(:,ttt) = iter_dat./dm;
  end


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

function [K] = gen_kernel(Tau,T,IRF_flag,resolution)

 Kernel = @(time,TimeConst)(((1./TimeConst).*(exp(- time * (1./ TimeConst)))));
 
    switch IRF_flag
        
        case 0 
            C1 = Kernel (Tau,T);
                K = C1;         
        case 1
            resolution = resolution;
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

    %   make Kernel and SVD and convolute

        C1 = Kernel (Tau,T);

            TK1 = zeros(ndata,L_basis);

            for con = 1:L_basis
                temp_conv1 = conv(C1(:,con),sim_irf);
                TK1(:,con) = temp_conv1(1:ndata);

            end
                K = TK1;

        
    end
end


function kld = kl(f,g)
        lr = log(f./g);
        lr(isnan(lr)|isinf(lr))=0;
        
        kld = (sum(f.*lr));      
end

function ir = irad(f,g)
        ir = kl(f,(f+g)/2) + kl(g,(f+g)/2);
end
