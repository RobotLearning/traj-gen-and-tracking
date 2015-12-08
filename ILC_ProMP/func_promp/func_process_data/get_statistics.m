function [data] = get_statistics(data)



    nSize = size(data(1).q, 2);
    nStates = size(data, 2);

    for j=1:nStates
        data(j).q_mean = mean(data(j).q);
        data(j).q_cov  = cov(data(j).q);
        data(j).q_var  = diag(data(j).q_cov);
        
        data(j).qdot_mean = mean(data(j).qdot);
        data(j).qdot_cov  = cov(data(j).qdot);
        data(j).qdot_var  = diag(data(j).qdot_cov);
        
        data(j).nSize  = nSize;
    end
    
    
    
end