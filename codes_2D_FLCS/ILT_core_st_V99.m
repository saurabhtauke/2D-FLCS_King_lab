function [F_ILT, parameter, Fitdata,NoiseStd,oppar] = ILT_core_st_V99(inputData,inAlpha,U1,S1,V1,U2,S2,V2,cond_num)
tic
    theMaxPoint = max(max(inputData));
    inData = inputData ./ theMaxPoint;
%         inData = inputData;

    if (inAlpha < 0 && inAlpha ~= -1 && inAlpha ~=-2 && inAlpha ~=-4)  
        disp 'Input alpha is not valid.'
        return
    end
    
    switch inAlpha
        case 0                
           disp 'BRD'
            Alpha_Auto = 1;
            AlphaStart = 1E3;
        case -4
            disp 'cascade'
            Alpha_Auto = 4;
            AlphaStart = 1E3;
        otherwise
            disp 'constant alpha'
            Alpha_Auto = 0;    
            AlphaStart = inAlpha;
    end
 
     ktemp = U1*S1*V1';
    ConditionNumber = cond_num;
    S = diag(S1)./ S1(1,1);
    x = find(S*ConditionNumber>1);
 	U1 = U1(:,x);
    S1 = S1(x,x);
    V1 = V1(:,x);

    S = diag(S2) ./ S2(1,1);
    x = find(S*ConditionNumber>1);
 	U2 = U2(:,x);
    S2 = S2(x,x);
    V2 = V2(:,x);

    K1 = U1 * S1 * V1';
    K2 = U2 * S2 * V2';
    
     % compressed data
    
     mtilde = U1' * inData' * U2;
    
    % projected data in range space
    mtilde2 = U1 * mtilde * U2';
    
%     NoiseStd = 1;
     NoiseStd = mean(std(((mtilde2 - inData'))))
%        NoiseStd = mean(std((sqrt(mtilde2) - sqrt(inData'))))
%        NoiseStd = (3/8)*mean(std(((mtilde2 - inData')./sqrt(mtilde2))))
%         NoiseStd = mean(std(( (mtilde2 - inData')./sqrt(inData'))))

        
        nn = (((mtilde2 - inData'))./sqrt(mtilde2));
         nn(isinf(nn) | isnan(nn)) = 0;

    w = (1./(mtilde2*theMaxPoint));
    w(isinf(w)) = 0;
    wvar = (numel(mtilde2))/sum(sum(w));

 
    switch Alpha_Auto 
        
        case {0,1,3,4}	

            [F_ILT, CompressedData, Chi, Alpha,oppar] = ...
                C_solver(inData,  U1, U2, V1, V2,...
                    S1, S2,AlphaStart, NoiseStd, Alpha_Auto, ConditionNumber,wvar,theMaxPoint);

        otherwise 
            	'invalid alpha'
    end
   oppar.time =  toc;
     
        F_ILT = F_ILT*theMaxPoint;
        Chi = Chi*theMaxPoint;

        parameter.chi = Chi;
		parameter.alpha = Alpha;
        Fitdata = K2*F_ILT*K1';
        oppar.comp_dat = CompressedData';
        oppar.NoiseStd = NoiseStd;
        oppar.mtilde = mtilde;
        oppar.mtilde2 = mtilde2;
        oppar.theMaxPoint = theMaxPoint;
        
end