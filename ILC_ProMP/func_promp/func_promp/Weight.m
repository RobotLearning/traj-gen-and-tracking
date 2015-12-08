
classdef Weight < handle

      
   properties
      nBasis;
      nJoints;
      nTraj;
      nDemo;
      
      demoq;
      w;
      mean; 
      cov;
      covchol;
      my_linRegRidgeFactor;
      
      w_full;
      mean_full;
      cov_full;
      covii_full;
      covchol_full;
      
      index;
   end
   
        
    
methods

    function obj = Weight(nBasis, nJoints, nTraj, nDemo)
        obj.nBasis = nBasis;
        obj.nJoints = nJoints;
        obj.nTraj = nTraj;
        obj.nDemo = nDemo;
    end

    function obj = least_square_on_weights(obj, Gn, demoq)
        % This computes the weights for the trajectory
        
        obj.demoq = demoq;
        dbg       = 0;
        
        nDemo    = obj.nDemo;
        nJoints  = obj.nJoints;
        nBasis   = obj.nBasis;

        my_linRegRidgeFactor = obj.my_linRegRidgeFactor ;
        
        % find weights by normal solution with Moore-Penrose Pseudo-Inverse
        MPPI = (Gn'*Gn  + my_linRegRidgeFactor) \ Gn';      
        
        % First get the weights of each trajectory for each joint
        for j=1:nJoints
            w{j}=[];
            for k = 1:nDemo
                w_   = MPPI*demoq{k}(:,j);
                w{j} = [w{j} ; w_'];

                if dbg % reproduce each trajectory
                    figurew(['learned weights ', num2str(k) , 'joint', num2str(j)]);
                    title 'Comparing trajectory reproduction by weights'
                    plot(demoq{k}(:,j), SBLUEBALLW(10));
                    plot(w_'*Gn', SREDBALL(5));
                end
            end
            ind{j} = [(j-1)*nBasis+1:1:j*nBasis];
        end
        w = w';
        obj.index = ind;

        % ============================================================
        % This is the important trick for the interaction!!!
        % We put all weights of all joints in a single matrix ...
        for j=1:nJoints
            obj.w_full = [obj.w_full w{j}];
        end
      
        % and compute the covariance matrix of all joints together!
        % This will link all the joints, from the human and from the robot
        % by means of their covariances
        obj.cov_full   = cov(obj.w_full);
        
        
        obj.mean_full  = mean(obj.w_full)';        
        
        % improve covariance
        [obj.cov_full, obj.covchol_full] = ...
             post_processing_covariance(obj.cov_full, nJoints*nBasis, nDemo);
         
    end % learn weights
        
end %methods


end






