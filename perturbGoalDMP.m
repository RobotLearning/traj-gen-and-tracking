% Randomly perturbs reference DMPs to a new goal position 

function perturbGoalDMP(dmp,sigma)

    dim = length(dmp);
    e = sigma * randn(dim,1);

    for i = 1:dim
        % perturb initial position slightly
        g = dmp(i).goal + e(i);
        dmp(i).setGoal(g);
    end            
    
        
end