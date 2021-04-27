function  [F_ILT, CompressedData, Chi, alpha,oppar] = C_solver(Data, ...
U1, U2, V1, V2,S1, S2, alpha, NoiseStd, flag, CondNum,wvar,theMaxPoint)

kt1 = U1*S1*V1';
kt2 = U2*S2*V2';

Number_T1 =length(V1); 
Number_T2 =length(V2);

n = 1;

  V = zeros((Number_T1)*Number_T2,1);
  
    for i = 1:length(S2)
        for j = 1:length(S1)
            if (S2(i,i)*S1(j,j) > S1(1,1)*S2(1,1)/CondNum)
				S(n)= S2(i,i) * S1(j,j);
          	CompressedData(n) = U2(:,i)'*Data*U1(:,j);
           CompressedNoiseVar(n) = 1;
          	V(:,n) = kron(V2(:,i), V1(:,j));
          	n = n+1;
            end
        end
    end

%  e= 1;  
 S = diag(S);
 n = length(S);
 data = CompressedData';
 
 fprintf(1, 'Size of compressed data = %d\n',n);
  
% modified K
K = S*V';
oppar.kk = K;

% optim options
options = optimset; options.GradObj = 'on'; options.Hessian = 'on';
options.LargeScale = 'on'; options.MaxIter = 5000;
options.Display = 'final';
options.TolX = 1e-8;
options.TolFun = 1e-16;
options.MaxFunEvals = 5000;


Identity = eye(n, n);

C_hat = ones(n, 1);

st_alpha = alpha;
init_alpha = 1E6;
pp=1;
alpha = init_alpha;
dec_rate = 0.2; %0.9 for smooth


stopper = 1/sqrt(NoiseStd)
oppar.stopper= stopper;
bias = 1000;
fstop = 0; 

while fstop ==0
	fprintf(1, ' %2.2e ', alpha);	

            [C_hat, fval, exitflag, output,grd] = fminunc('C_minunc', C_hat, options, ...
	    data,  K, alpha,Identity) ;

    fronorm =  norm(C_hat, 'fro');
	Chi = alpha *fronorm/stopper;
    oppar.fronorm(pp) = fronorm;
    oppar.chi(pp) = Chi;
    oppar.fval(pp) = -fval;
    
        switch flag
            case 0  % fixed
                
                alpha = alpha * dec_rate;
                oppar.cv_dump(:,pp) = C_hat;
                oppar.apar(pp) = alpha;
                
                pp = pp+1;
                
                if(fronorm  > stopper) 
%                     fstop =1; 
                    'chi_stop'
                end
                
                
                if  alpha<st_alpha
                    fstop = 1;
                    'min alpha reached'
                end
            case 1  % BRD
%                 
                    alpha_opt =  sqrt(n).*(NoiseStd)./(norm(C_hat,'fro'));
                  		
                    oppar.cv_dump(:,pp) = C_hat;
                    oppar.apar(pp) = alpha_opt;
                    pp = pp+1;

                if (abs((alpha-alpha_opt)/alpha) < 1E-3) 
                    fstop = 1; 
                    'alpha_converged'
                end
                
                if(fronorm  > stopper*bias) 
                    fstop =1; 
                    'chi_stop'
                end
                alpha = alpha_opt;
               
              
                
            case 4
                
                alpha = alpha * dec_rate;
                oppar.cv_dump(:,pp) = C_hat; 
                pp = pp+1;
                
                if pp > 20 || alpha<1E-7
                    fstop = 1;
                end
                
            otherwise
                fstop = 1;
        end
        
        % Convert vector back into matrix
	F_ILT = max(0, reshape(K'*C_hat, (Number_T1), Number_T2)');	
    
end	

end	

