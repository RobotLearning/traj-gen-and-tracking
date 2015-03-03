% 
% Iterative Learning Control 
%
% updates trajectories using feedback
%
% includes ancillary methods that update 
% i) trajectory directly
% ii) weights of the trajectory
% iii) weights of the derivative of the trajectory
%
% switching between different methods done in initialization
%

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
        
        % update method
        upd_meth
        % last input - trj, or weights
        inp_last
        % output matrix
        C
        % Psi matrix used for weight updates
        Psi
        % learning matrix
        L
    end
    
    methods
        
        function obj = tILC(traj,model,update_method)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'tILC';
            obj.error = 0;
            obj.upd_meth = update_method;
            obj.formL(traj);
            
            N = traj.N - 1;
            h = traj.t(2) - traj.t(1);
            obj.C = model.C;            
            Cout = obj.C;
            dim_x = size(Cout,2);
            obj.Psi = obj.formPsi(dim_x*N,traj.t);
            
            sbar = Cout'*((Cout*Cout')\traj.s);
            svec = sbar(:);

            switch update_method 
                case 't'                
                obj.inp_last = svec(1:end-dim_x);
                case 'w' 
                svec = sbar(:);
                obj.inp_last = pinv(obj.Psi) * svec(1:end-dim_x);
                case 'wd'
                svec = sbar(:);
                D = eye(N*dim_x) - diag(ones(1,(N-1)*dim_x),-dim_x);
                obj.inp_last = pinv(obj.Psi) * D * svec(1:end-dim_x);
                otherwise 
                error('Something wrong!');
            end
            
            
        end
        
        % Form the Psi matrix
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
        
        % Form the learning matrix
        function formL(obj,traj)
            
            N = traj.N - 1;
            K = traj.K;
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
            
            obj.L = pinv(Kl)*L;
            
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
            Cout = obj.C;
            dim = size(Cout,2);
            dev = y - traj.s;   
            
            % ILC update
            w = obj.inp_last + pinv(obj.Psi)*obj.L*dev(:);
            obj.inp_last = w;
            % form back the trajectory
            sbar = obj.Psi * w;
            sbar = reshape(sbar,dim,N);
            sbar = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t, Cout*sbar, traj.unom, K);
            
        end
        
        % doing regression on derivatives is equivalent
        function traj2 = updateDerivativeWeights(obj,traj,y)
            
            h = traj.t(2) - traj.t(1);
            N = traj.N - 1;
            K = traj.K;
            dev = y - traj.s; 
            Cout = obj.C;
            
            % construct derivative matrix
            dim = size(Cout,2);
            D = eye(N*dim) - diag(ones(1,(N-1)*dim),-dim);
            
            % construct integral matrix S
            Sv = repmat([1;0],N,1);
            Sv2 = repmat([0;1],N,1);
            S = repmat([Sv,Sv2],1,N);
            S = tril(S);
            
            % ILC update
            wd = obj.inp_last + pinv(obj.Psi)*D*obj.L*dev(:);
            obj.inp_last = wd;
            % form back the trajectory
            sbar = S * obj.Psi * wd;
            sbar = reshape(sbar,dim,N);
            sbar = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t, Cout*sbar, traj.unom, K);
            
            % get s from derivatives
            %S = h*tril(ones(N+1));
            %sbar = S * sd;           

            
        end
        
        % update trajectories 
        function traj2 = updateTraj(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - traj.s;
            Cout = obj.C;
            dim = size(Cout,2);      

            % ILC update
            s = obj.inp_last + obj.L*dev(:);
            obj.inp_last = s;
            % form back the trajectory
            sbar = reshape(s,dim,N);
            sbar = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t, Cout*sbar,traj.unom,K);
            
        end
        
        % switching function that switches between different
        % methods
        function traj2 = feedforward(obj,traj,y)
        
            update_method = obj.upd_meth;
            
            switch update_method 
                case 't'
                traj2 = updateTraj(obj,traj,y);
                case 'w'
                traj2 = updateWeights(obj,traj,y);
                case 'wd'
                traj2 = updateDerivativeWeights(obj,traj,y);
                otherwise %TODO
                    %TODO
            end
            
        end
    end
    
end