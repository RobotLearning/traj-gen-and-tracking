function [visited, joint_KF] = mainloop(pmp, prm, obs)

    h_temp_obs   = [];   visited.k    = [];     visited.q    = [];
    visited.qdot = [];   visited.mse  = [];
    
    sigma_obs        = obs.stdev;
    
    % feature: Gaussian
    P0   = pmp.w.cov_full;
    x0   = pmp.w.mean_full;
    
    % The observation noise is nJoints*2 as it also includes the velocity
    Robs_pmp  = (sigma_obs^2)*eye( 2*(prm.nJointA+prm.nJointB));
    
    % Recursive processing =====================
    aux.nObsTotal = length(obs.index);
    
    dtPlot = round(aux.nObsTotal/prm.plot.nPlots);
    plot_at_k = 1:dtPlot:aux.nObsTotal;
    if plot_at_k(end)~=aux.nObsTotal
        plot_at_k(end+1)=aux.nObsTotal;
    end
    
    
    aux.counter = 1;
    
    for k = obs.index
        
        visited.k = [visited.k                k];
        visited.q = [visited.q  obs.q(:,k)];
        visited.qdot = [visited.qdot  obs.qdot(:,k)];        
  

        try delete(pmp.h.qshade);   pmp.h.qshade =[]; end

        % create observation matrix
        H0   = observation_matrix(k, pmp, obs.joint, obs.jointvel);        
        
        z0 = [];
        for j=1:pmp.nJoints
            q_qdot = [obs.q(j, k); obs.qdot(j, k)];
            z0  = [z0; q_qdot];
        end
       
        % Gaussian conditioning using Kalman filtering
        % Equation (5) and (6) of the ProMP paper (NIPS2013)
        [x1, P1] = KF_recursion(x0, P0, H0, z0, Robs_pmp);
        
        regulariz_demo = pmp.nDemo;
        % regulariz_demo = 100;
        [P1, ~] = post_processing_covariance(P1, pmp.nJoints*pmp.nBasis, regulariz_demo);
        x0 = x1;
        P0 = P1;

        % separate the huge full matrices per each joint (most
        % useful to plot, but useless for the KF itself)
        joint_KF = perJoint_pmp(x1, P1, pmp);
        
        % compute current prediction error
        for j=1:prm.nJointB
            error_u(j) = RMSE(pmp.t, [joint_KF(prm.nJointA+j).q_mean  joint_KF(prm.nJointA+j).q_mean(end)],...
                 pmp.t, obs.groundtr(prm.nJointA+j).q);
        end
        visited.mse = [visited.mse  mean(error_u)];
        
        
        fprintf('%g  observations left\n', aux.nObsTotal-aux.counter+1);
        
        
        if 0% sum(aux.counter == plot_at_k)
            
            pmp.h = plotTypeIterative(pmp.h, pmp.t, pmp.nTraj, joint_KF, prm);
            pmp.h = plot_observations(pmp.h, pmp.t, prm.nJointA, obs, visited);
            
            if  0%prm.saveSVG.at_every_conditioniong_plot
                plot2svg(['task' PMP{p}.pair '_A_iter_' num2str(aux.counter) '.svg']);
            end                          

            drawnow;
        end
        aux.counter = aux.counter+1;
    end
                                      
    if  1% sum(aux.counter == plot_at_k)

        pmp.h = plot_observations(pmp.h, pmp.t, prm.nJointA, obs, visited);        
        pmp.h = plotTypeIterative(pmp.h, pmp.t, pmp.nTraj, joint_KF, prm);

                  

        drawnow;
    end
                     

end


