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
            
            N = traj.N - 1;
            h = traj.t(2) - traj.t(1);
            obj.trj = traj.s;
            %s = dmp.evolve();
            %s = s(1,:);
            C = model.C;
            obj.sbar = C'*((C*C')\traj.s);
            dim = size(obj.sbar,1);
            
            obj.Psi = obj.formPsi(dim*N,traj.t);
            sbar = obj.sbar(:);

            obj.w = pinv(obj.Psi) * sbar(1:end-dim);

            D = eye(N*dim) - diag(ones(1,(N-1)*dim),-dim);
            obj.wd = pinv(obj.Psi) * D * sbar(1:end-dim);
            
        end
        
        function Psi = formPsi(obj,bfs,t)
            N = length(t)-1;
            Psi = zeros(2*N,bfs);
            h = ones(bfs,1) * bfs^(1.5);
            c = linspace(t(1),t(end),bfs);
            tbar = linspace(t(1),t(end),2*N);
            for i = 1:bfs
                Psi(:,i) = obj.basis(tbar,h(i),c(i));
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
    
            % get rid of x0 in dev
            %ddev = diff(dev')';
            %dev = dev(:,2:end);
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            %delta = a_p * dev + a_d * ddev;
            % construct learning matrix
            L = diag((a_p+a_d)*ones(1,N),1) - a_d*eye(N+1);
            L = L(1:end-1,:);            
            
            Kl = zeros(size(K,1)*N,size(K,2)*N);
            % construct lifted domain matrix Kl
            for l = 1:N
                vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                Kl(vec1,vec2) = K(:,:,l); % on diagonals only
            end
            
            obj.w = obj.w + pinv(obj.Psi)*pinv(Kl)*L*dev(:);
            sbar = obj.Psi * obj.w;
            dim = size(obj.sbar,1);
            sbar = reshape(sbar,dim,N);
            obj.sbar = [sbar,obj.sbar(:,end)];
            
            traj2 = Trajectory(traj.t, obj.sbar, traj.unom, K);
            
        end
        
        % TODO: show that doing regression on derivatives is equivalent
        function traj2 = updateDerivativeWeights(obj,traj,y)
            
            h = traj.t(2) - traj.t(1);
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj; 
            
            % construct derivative matrix
            dim = size(obj.sbar,1);
            D = eye(N*dim) - diag(ones(1,(N-1)*dim),-dim);
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            %delta = a_p * dev + a_d * ddev;
            % construct learning matrix
            L = diag((a_p+a_d)*ones(1,N),1) - a_d*eye(N+1);
            L = L(1:end-1,:);            
            
            Kl = zeros(size(K,1)*N,size(K,2)*N);
            % construct lifted domain matrix Kl
            for l = 1:N
                vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                Kl(vec1,vec2) = K(:,:,l); % on diagonals only
            end
            
            % construct integral matrix S
            Sv = repmat([1;0],N,1);
            Sv2 = repmat([0;1],N,1);
            S = repmat([Sv,Sv2],1,N);
            S = tril(S);
            
            obj.wd = obj.wd + pinv(obj.Psi)*D*pinv(Kl)*L*dev(:);
            sbar = S * obj.Psi * obj.wd;
            sbar = reshape(sbar,dim,N);
            obj.sbar = [sbar,obj.sbar(:,end)];
            
            traj2 = Trajectory(traj.t, obj.sbar, traj.unom, K);
            
            % get s from derivatives
            %S = h*tril(ones(N+1));
            %sbar = S * sd;           

            
        end
        
        % update trajectories 
        function traj2 = feedforward(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - obj.trj;
            %dev = dev(2,:);
            %h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.2;
            
            % construct learning matrix
            L = diag((a_p+a_d)*ones(1,N),1) - a_d*eye(N+1);
            L = L(1:end-1,:);            
            
            Kl = zeros(size(K,1)*N,size(K,2)*N);
            % construct lifted domain matrix Kl
            for l = 1:N
                vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                Kl(vec1,vec2) = K(:,:,l); % on diagonals only
            end
                        
            dim = size(obj.sbar,1);
            sbar = obj.sbar(:);
            s = sbar(1:end-dim);
            s = s + pinv(Kl)*L*dev(:);
            sbar = reshape(s,dim,N);
            obj.sbar = [sbar,obj.sbar(:,end)];
            
            traj2 = Trajectory(traj.t, obj.sbar,traj.unom,K);
            
        end
        
    end
    
end