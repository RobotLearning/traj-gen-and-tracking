function [H] = observation_matrix(k, p, observedJointPos, observedJointVel)
% Generic implementation


    nJoints = p.nJoints;
    nTraj   = p.nTraj; % number of steps in the trajectory for EACH joint
    
    Gn   = p.basis.Gn;
    Gn_d = p.basis.Gndot;
                 
    % introduce normalizedTime
    % normalizedTime = round( timepoint / dt );
    normalizedTime = k;
 
    % return the value of all bases at the specified normalizedTime
    Hq_measured    = Gn(normalizedTime,:); 
    Hqdot_measured = Gn_d(normalizedTime,:);
    
    % for joints that are not being measured the basis are zero
    Hq_unmeasured    = zeros(1, p.nBasis);
    Hqdot_unmeasured = zeros(1, p.nBasis);
    
    
    % create the full observation matrix
    % fprintf('***Creating obs. matrix\n')
    for j=1:nJoints
        %fprintf('\nJoint %g\n', j);
        
        % adding features for position
        if  isempty(observedJointPos)
            %fprintf('No Position features\n')
            H_temp = Hq_unmeasured;
        else
            if sum(j==observedJointPos)==0
                %fprintf('No Position features\n')
                H_temp = Hq_unmeasured;
            else
                %fprintf('Position features added\n')
                H_temp = Hq_measured;                
            end
        end
        
        % adding features for velocity
        if isempty(observedJointVel)
            %fprintf('No velocity features\n')
            H_temp = [H_temp; Hqdot_unmeasured];
        else
            if sum(j==observedJointVel)==0
                %fprintf('No velocity features\n')
                H_temp = [H_temp; Hqdot_unmeasured];
            else
                %fprintf('Velocity features added\n')
                H_temp = [H_temp; Hqdot_measured];                
            end  
        end
        H{j} = H_temp;
    end

    
    % create the block diagonal matrix
    H = blkdiag( H{:} );
    
    
    
    
end
    
    