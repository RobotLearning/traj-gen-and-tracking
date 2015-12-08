function [] = animate_bioRob(bioRob, q, nTraj)



    [m,n] = size(q);
    
    if m >= nTraj
        
        t_orig  = linspace(0,1,m);
        t_small = linspace(0,1,nTraj);
        q_small = interp1(t_orig, q, t_small);
        
    else
        q_small = q;
    end
    
     bioRob.plot(q_small);

end