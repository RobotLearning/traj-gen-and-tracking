% Randomly perturbs reference DMPs to a new initial condition (position and
% velocity)

function [Q0,E0] = perturbInitCondDMP(dmp,sigma)

    dim = length(dmp);
    e0 = sigma * randn(dim,1);
    ed0 = sigma * randn(dim,1);

    for i = 1:dim
        % perturb initial position slightly
        q0(i) = dmp(i).y0(1) + e0(i);
        qd0(i) = dmp(i).y0(2) + ed0(i);
        % append zero acc
        y0 = [q0(i);qd0(i);0];
        dmp(i).y0 = y0;
        dmp(i).resetStates();
    end            
    
    Q0 = [q0(:);qd0(:)];
    E0 = [e0(:);ed0(:)];
        
end