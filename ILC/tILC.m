% (basic) Iterative Learning Control using PD-type trajectory update

classdef tILC < ILC
    
    % fields common to all ILCs (defined in abstract ILC class)
    properties
         
        % number of total episodes so far
        episode
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        error
        
        % ILC's Last input sequence (in this case a trajectory)
        u_last
        % initial trajectory's dimension increased
        sbar
    end
    
    methods
        
        function obj = tILC(traj,model)
                        
            obj.episode = 0;
            obj.color = 'b';
            obj.name = 'ILC modifying trajectories';
            obj.error = 0;
            
            N = traj.N - 1;
            obj.u_last = traj.s;
            %s = dmp.evolve();
            %s = s(1,:);
            C = model.C;
            obj.sbar = C'*((C*C')\traj.s);
            
        end
        
        % update weights for the trajectory instead
        function traj2 = updateWeights(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.u_last;
            %dev = dev(2,:);
            %h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            ddev = diff(dev')';
            dev = dev(:,2:end);
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            delta = a_p * dev + a_d * ddev;
            
            for i = 1:N
                obj.sbar(:,i) = obj.sbar(:,i) + pinv(K(:,:,i)) * delta(:,i);
            end
            
            traj2 = Trajectory(traj.t, obj.sbar,traj.unom,K);
            
        end
        
        % update trajectories 
        function traj2 = feedforward(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.u_last;
            %dev = dev(2,:);
            %h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            ddev = diff(dev')';
            dev = dev(:,2:end);
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            delta = a_p * dev + a_d * ddev;
            
            for i = 1:N
                obj.sbar(:,i) = obj.sbar(:,i) + pinv(K(:,:,i)) * delta(:,i);
            end
            
            traj2 = Trajectory(traj.t, obj.sbar,traj.unom,K);
            
        end
        
        % override ilc's record method
        function record(obj,u_applied,cost)            
            obj.episode = obj.episode + 1;
            obj.u_last = u_applied;
            obj.error(obj.episode) = cost;            
        end
        
    end
    
end