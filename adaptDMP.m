% Adapts reference DMPs to a new goal position
function [dmpNew,s] = adaptDMP(yin,goal,dmpOrigin,wOrigin,steps)

    dim = length(dmpOrigin);
    can = dmpOrigin(1).can;
    alpha = dmpOrigin(1).alpha_g;
    beta = dmpOrigin(1).beta_g;
    numbf = size(wOrigin,2);
    

    for i = 1:dim
        % append zero velocity
        y0 = [yin(i);0;0];
        % create the dmp trajectory
        dmpNew(i) = DDMP(can,alpha,beta,goal(i),y0,numbf);
        % set the original weights
        dmpNew(i).w = wOrigin(i,:);
        % evolve the DMP
        [x,si] = dmpNew(i).evolve(steps);         
        s(i,:) = si(1,:);

    end 
end