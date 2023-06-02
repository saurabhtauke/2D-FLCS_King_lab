function [map_ab] = partition_2D_core(n_gates,del_lb,del_ub,a_macro_time,a_gate,b_macro_time,b_gate)
    
     map_ab = zeros(n_gates,n_gates);
  
    for i = 1 : length(a_macro_time)
       sample_T = a_macro_time(i);
       sample_gate = a_gate(i);
       temp_ab  = ((b_macro_time>sample_T+del_lb)&(b_macro_time<sample_T+del_ub));  
       map_y12 = b_gate(temp_ab);
             
       if(~isempty(map_y12))
            for ii = 1:length(map_y12)
                map_ab(sample_gate,map_y12(ii)) = map_ab(sample_gate,map_y12(ii))+1;
            end
       end       
    end
end