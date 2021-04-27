function [f_val,grd,Hss] = C_minunc(c, data,  K, alpha,Identity)  

% V3 done 
% before G
    
        ktc = K'*c;
        mx = max(0,ktc);
        
        f_val = 0.5*((ktc'*mx)+ c'*c*alpha) - c'*data;
        
        grd = ((K*mx)+ alpha.*c) - data;
        
        Hss = (K*max(0,K')+alpha);
   

end
