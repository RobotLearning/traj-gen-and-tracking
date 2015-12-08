function u = torque_controller(robot, t, q, qdot, varargin)
%
%
%
%
% 24.04.2014: Implenentation of the controller function. v00. Guilherme
    

    opt        = varargin{1};
    
    %% Compute the error 
    % find the reference value for a given time t (t is controlled by
    % ode45)
    q_ref_t = interp1(opt.ref.t, opt.ref.q, t);
    qerror     = q_ref_t - q;
    
    % qerrordot  = opt.ref.qdot - qdot; Velocity controller is only damping
    qdot_ref_t = interp1(opt.ref.t, opt.ref.qdot, t);
    qerrordot    =  qdot_ref_t - qdot; 
      
    
    %% Add control commands to u
    u = 0;
    
    if opt.controller.feedback
        
        % compute pd commands for each joint of the arm
        u   =     qerror.*opt.feedback.P;
        u   = u + qerrordot.*opt.feedback.D;
        
        % integral gain implementation not ready yet but can be implemented
    end
    
    if opt.controller.feedforward
        % feedforward commands are given by ILC input
        u = u + interp1(opt.ref.t, opt.ilc.u, t);
    end

end