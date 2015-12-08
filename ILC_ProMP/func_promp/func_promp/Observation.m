classdef Observation < handle

      
   properties

      joint;
      jointvel;
      index;
      stdev;
      q;
      qdot;
      groundtr; % observed ground truth
      groundtr_names;
      
      nTotalJoints;
      nTraj;
      
   end
   
        
    
methods

    function obj = Observation(stdev, param, groundtr)

        obj.joint    = param.obs.jointPos;
        obj.jointvel = param.obs.jointVel;
        obj.stdev    = stdev;
        
        
        obj.nTotalJoints  = size(groundtr,2);
        obj.nTraj = numel(groundtr(1).q);
        
        if 0
            temp=[];
            for j=1:obj.nTotalJoints
                temp = [temp   groundtr(j).q'     groundtr(j).qdot'];
            end
            obj.groundtr = temp;
            obj.groundtr_names ={'q1 q1dot q2 q2dot q3 q3dot u1 u1dot u2 u2dot u3 u3dot'};            
        else
            obj.groundtr = groundtr;
        end
        

        
        
        
        % prepare the vector of observations by filling it with zeros
        for j = 1:obj.nTotalJoints
            obj.q(j,:)    = 0.*ones(1, obj.nTraj);
            obj.qdot(j,:) = 0.*ones(1, obj.nTraj);
        end
        
        
        
    end

    function obj = measuredIndexes(obj, index)
        
        
        obj.index = index;

        if max(obj.index)>obj.nTraj
            error('Observing more points than points available in the test data');
        end

        for j = obj.joint
            obj.q(j,index)    = obj.groundtr(j).q(index);            
        end
        
        for j = obj.jointvel
            obj.qdot(j,index) = obj.groundtr(j).qdot(index); 
        end
        
        
% % %         for j = obj.joint
% % %             obj.q(j,index)    = obj.groundtr(index,j);            
% % %         end
% % %         
% % %         % velocity is usually not observed but I put it here just for
% % %         % completeness
% % % % % %         qdot = diff(obj.groundtr);
% % % % % %         qdot = [qdot; qdot(end,:)];
% % %         
% % %         for j = obj.jointvel
% % %             ii = 2*(j-1)+2;
% % %             obj.qdot(j,index) = obj.groundtr(index, ii );
% % %         end

    end

    
end %methods


end


%     clear obs
%     obs(1).parent.set    = 2;
%     obs(1).parent.joints = [1 2];
%     obs(1).parent.time   = [   [50:55]   [65]   80   90    [10:15]];
%     obs(1).parent.stdev  = 0.01;
%     for j = obs(1).parent.joints
%         index_observed_data = obs(1).parent.time;
%         data_set = obs(1).parent.set;
%         obs(j).q = NaN.*ones(1,data(1).nSize) ;
%         obs(j).q(index_observed_data ) = data(j).q(data_set, index_observed_data );
%     end
