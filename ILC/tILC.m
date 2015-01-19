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
        % weights of the derivatives
        wd
    end
    
    methods
        
        function obj = tILC(traj,model)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'tILC';
            obj.error = 0;
            
            N = traj.N;
            h = traj.t(2) - traj.t(1);
            obj.trj = traj.s;
            %s = dmp.evolve();
            %s = s(1,:);
            C = model.C;
            obj.sbar = C'*((C*C')\traj.s);
            
            obj.Psi = obj.formPsi(N,traj.t);
            obj.w = pinv(obj.Psi) * obj.sbar';

            D = (diag(-ones(1,N-1),-1) + eye(N))/h;
            obj.wd = pinv(obj.Psi) * D * obj.sbar';
            
        end
        
        function Psi = formPsi(obj,bfs,t)
            N = length(t);
            Psi = zeros(N,bfs);
            h = ones(bfs,1) * bfs^(1.5);
            c = linspace(t(1),t(end),bfs);
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
            
            h = traj.t(2) - traj.t(1);
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj;            
    
            % construct learning matrix
            L = 0.5;
            % get rid of x0 in dev
            ddev = diff(dev')';
            dev = dev(:,2:end);
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            delta = a_p * dev + a_d * ddev;
            
            Kl = zeros(size(K,1)*N,size(K,2)*N);
            % construct lifted domain matrix Kl
            for l = 1:N
                vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                Kl(vec1,vec2) = K(:,:,l); % on diagonals only
            end
            
            obj.w = obj.w + pinv(obj.Psi)*pinv(Kl)*L*dev(:);
            sbar = obj.Psi * obj.w;
            
            traj2 = Trajectory(traj.t, sbar', traj.unom, K);
            
        end
        
        % TODO: show that doing regression on derivatives is equivalent
        function traj2 = updateDerivativeWeights(obj,traj,y)
            
            h = traj.t(2) - traj.t(1);
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj; 
            
            % set learning rate
            a = 0.5;
            
            Kl = zeros(size(K,1)*N,size(K,2)*N);
            % construct lifted domain matrix Kl
            for l = 1:N
                vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                Kl(vec1,vec2) = K(:,:,l); % on diagonals only
            end
            
            % get derivatives
            D = (diag(-ones(1,N),-1) + eye(N+1))/h;
            S = h*tril(ones(N+1));
            obj.wd = obj.wd + a*pinv(obj.Psi)*D*dev(:);
            sd = obj.Psi * obj.wd;
            sbar = S * sd;           
            
            traj2 = Trajectory(traj.t, sbar', traj.unom,K);
            
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