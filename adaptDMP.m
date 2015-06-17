% Adapts reference DMPs to a new goal position
function [dmpNew,s] = adaptDMP(yin,goal,dmp,steps)

    dim = length(dmp);
    can = dmp(1).can;
    alpha = dmp(1).alpha_g;
    beta = dmp(1).beta_g;
    numbf = size(dmp.w,2);
    

    for i = 1:dim
        % append zero velocity
        y0 = [yin(i);0;0];
        % create the dmp trajectory
        dmpNew(i) = DDMP(can,alpha,beta,goal(i),y0,numbf);
        % set the original weights
        dmpNew(i).w = dmp.w(i,:);
        % evolve the DMP
        [x,si] = dmpNew(i).evolve(steps);         
        s(i,:) = si(1,:);

    end 
end