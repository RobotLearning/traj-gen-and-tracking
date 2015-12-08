
function kf = perJoint_pmp(xFull, Pfull, pmp)
% This separates the full mean covariance matrix in the components of each
% joinit.

    
    nBasis = pmp.nBasis;
    nTraj = pmp.nTraj;
    
    for j=1:pmp.nJoints
        ind = [(j-1)*nBasis+1:1:j*nBasis];
        kf(j).w_mean        = xFull(ind);
        kf(j).w_Sigma       = Pfull(ind,ind);
        kf(j).w_Sigma_ii    = diag(  kf(j).w_Sigma  );
        
        % from weights to position and velocity
        [q_mean, q_Sigma_ii, qdot_mean, qdot_Sigma_ii] = ...
                   theta_to_traj(  kf(j).w_mean, kf(j).w_Sigma, ...
                                pmp.basis, pmp.phase.dt, nTraj );
        
        if 1
            kf(j).q_mean        = [q_mean];
            kf(j).q_Sigma_ii    = [q_Sigma_ii];

            kf(j).qdot_mean     = [qdot_mean];
            kf(j).qdot_Sigma_ii = [qdot_Sigma_ii];
            
        else  % add initial state to make size N and not N-1
            kf(j).q_mean        = [ q_mean];
            kf(j).q_Sigma_ii    = [  q_Sigma_ii];

            kf(j).qdot_mean     = [ qdot_mean];
            kf(j).qdot_Sigma_ii = [ qdot_Sigma_ii]; 
        end
                            
    end
    
end

function [x_mean, x_Sigma_ii, xdot_mean, xdot_Sigma_ii] = ...
                         theta_to_traj(w_mean, P_w, basis, phase_dt, nTraj)
% Given the weights and basis of the ProMP compute the trajectory.
         
    x_mean         = [];
    x_Sigma_ii     = [];
    xdot_mean      = [];
    xdot_Sigma_ii  = [];

    for i = 1:(nTraj-1)

        timePoint = i*phase_dt;
        [mu_x, ~, Sigma_t] =  get_distributions_at_time_T1(w_mean,...
                                          P_w, basis, phase_dt, timePoint);

        x_mean        = [x_mean          mu_x(1)];
        xdot_mean     = [xdot_mean       mu_x(2)];            
        x_Sigma_ii    = [x_Sigma_ii      Sigma_t(1,1)];        
        xdot_Sigma_ii = [xdot_Sigma_ii   Sigma_t(2,2)];

    end


end
    
    
    
    
    
    
    
    