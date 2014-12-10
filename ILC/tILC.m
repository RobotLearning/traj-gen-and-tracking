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
        
        % TODO: NOT USED!
        u_last
        % ILC's Last input sequence (in this case a trajectory)
        trj
        % initial trajectory's dimension increased
        sbar
        % Psi matrix used for weight updates
        Psi
        % weights instead of the trajectory
        w
    end
    
    methods
        
        function obj = tILC(traj,model)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'ILC modifying trajectories';
            obj.error = 0;
            
            N = traj.N - 1;
            obj.trj = traj.s;
            %s = dmp.evolve();
            %s = s(1,:);
            C = model.C;
            obj.sbar = C'*((C*C')\traj.s);
            
            obj.Psi = obj.formPsi(N,traj.t(1:N));
            s = traj.s(:,1:N);
            obj.w = pinv(obj.Psi) * s(:);
            
        end
        
        function Psi = formPsi(obj,bfs,t)
            N = length(t);
            Psi = zeros(N,bfs);
            h = ones(bfs,1) * bfs^(1.5);
            c = linspace(t(1),t(end-1),bfs);
            for i = 1:bfs
                Psi(:,i) = obj.basis(t,h(i),c(i));
            end
        end
        
        % basis functions are unscaled gaussians
        function out = basis(obj,t,h,c)
        out = exp(-h * (t - c).^2);
        end
        
        % update weights for the trajectory instead
        function traj2 = updateWeights(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj;
            dev = dev(:,2:end);
    
            % set learning rate
            a = 0.5;
            obj.w = obj.w + a*pinv(obj.Psi)*dev(:);
            s = obj.Psi * obj.w;
            s = [s; traj.s(:,end)]';
            
            traj2 = Trajectory(traj.t, s, traj.unom,K);
            
        end
        
        % update trajectories 
        function traj2 = feedforward(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj;
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
        
    end
    
end