% Differentiates DMPs

function [ddmp,s] = diffDMP(yin,goal,dmp)

    dim = length(dmp);
    can = dmp(1).can;
    alpha = dmp(1).alpha_g;
    beta = dmp(1).beta_g;
    numbf = length(dmp(1).FORCE.w);

    for i = 1:dim
        % append zero acceleration
        y0 = [yin(i);0];
        % create the dmp trajectory
        ddmp(i) = discreteDMP(can,alpha,beta,goal(i),y0,numbf);
        [~,ref] = dmp(i).evolve();
        dref = diff(ref(1,:))/can.dt;
        dref(end+1) = dref(end);
        % set the original weights
        ddmp(i).setWeights(dref);
        % evolve the DMP
        [x,si] = ddmp(i).evolve();         
        s(i,:) = si(1,:);

    end 
end