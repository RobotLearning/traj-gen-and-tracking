% 
% Iterative Learning Control 
%
% updates trajectories/DMPs using feedback
%
% includes ancillary methods that update 
% i) trajectory directly
% ii) weights of the trajectory
% iii) weights of the derivative of the trajectory
% iv) weights of the DMPs
%
% switching between different methods done in initialization
%

classdef wILC < ILC
    
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
        
        % flag for modifying only positions
        modifyPosOnly
    end
    
    methods
        
        %% Constructor for wILC
        function obj = wILC(traj,model,varargin)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'wILC';
            obj.error = 0;
            
            N = traj.N - 1;
            h = traj.t(2) - traj.t(1);
            obj.C = model.C;            
            Cout = obj.C;
            dim = size(Cout,2);
            obj.Psi = obj.formPsi(dim*N,traj.t);
            
            % modify only positions
            obj.modifyPosOnly = false;
            
            % learning matrix
            obj.formL(model,traj);
            
            sbar = Cout'*((Cout*Cout')\traj.s);
            svec = sbar(:);

            switch varargin{1}
                case 't'                
                obj.inp_last = svec(1:end-dim);
                obj.upd_meth = 't';
                case 'w' 
                obj.inp_last = obj.Psi \ svec(1:end-dim);
                obj.upd_meth = 'w';
                case 'wd'
                D = eye(N*dim) - diag(ones(1,(N-1)*dim),-dim);
                obj.inp_last = pinv(obj.Psi) * D * svec(1:end-dim);
                obj.upd_meth = 'wd';
                otherwise 
                %dmp = varargin{1};
                %dmp.constructF(traj.t);
                obj.inp_last = svec(1:end-dim);
                %obj.inp_last = dmp.FORCE.w;
                obj.upd_meth = 'dmp';
            end
            
            
        end
        
        %% Form the Psi matrix
        function Psi = formPsi(obj,bfs,t)
            N = length(t)-1;
            dimx = size(obj.C,2);
            Psi = zeros(dimx*N,bfs);
            h = ones(bfs,1) * bfs^(2.5);
            c = linspace(t(1),t(end),bfs);
            tbar = linspace(t(1),t(end),dimx*N);
            for i = 1:bfs
                Psi(:,i) = obj.basis(tbar,h(i),c(i));
            end
        end
        
        %% Form the learning matrix
        
        % Arimoto style PD-update
        function basicL(obj,traj)
            
            N = traj.N - 1;
            dim = size(obj.C,2)/2;
            dim = floor(dim);
            % set learning rate
            a_p = .5;
            a_d = .2;
            %delta = a_p * dev + a_d * ddev;
            % construct learning matrix
            %L = diag((a_p+a_d)*ones(1,N),1) - a_d*eye(N+1);
            
            M = [a_d*eye(dim),a_p*eye(dim)];
            M = kron(eye(N),M);
            obj.L = [zeros(dim*N,dim),M,zeros(dim*N,dim)];       
        end
        
        % Model based ILC update
        function lift(obj,model,traj)
            
            N = traj.N - 1;
            dim_x = model.SIM.dimx;
            dim_y = model.SIM.dimy;
            dim_u = model.SIM.dimu;
            F = zeros(N*dim_x, N*dim_u);
            
            % deal C matrix to G
            G = cell(1,N);
            Ql = cell(1,N);
            Rl = cell(1,N);
            % TODO: this could also be time varying!
            [G{:}] = deal(model.C);
            [Ql{:}] = deal(model.COST.Q);
            [Rl{:}] = deal(model.COST.R);
            G = blkdiag(G{:});
            Ql = blkdiag(Ql{:});
            Rl = blkdiag(Rl{:});
            
            if isa(model,'Linear')
            
                Ad = model.Ad;
                Bd = model.Bd;
                % construct lifted domain matrix F
                % TODO: this can be computed much more efficiently
                % for linear systems
                for i = 1:N
                    for j = 1:i
                        vec_x = (i-1)*dim_x + 1:i*dim_x;
                        vec_u = (j-1)*dim_u + 1:j*dim_u;
                        mat = Bd;
                        for k = j+1:i
                            mat = Ad * mat;
                        end
                        F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                    end
                end
            else
                % get linear time variant matrices around trajectory
                [Ad,Bd] = model.linearize(traj);
                
                % construct lifted domain matrix F
                for i = 1:N
                    for j = 1:i
                        vec_x = (i-1)*dim_x + 1:i*dim_x;
                        vec_u = (j-1)*dim_u + 1:j*dim_u;
                        mat = Bd(:,:,j);
                        for k = j+1:i
                            mat = Ad(:,:,k) * mat;
                        end
                        F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                    end
                end
            end
            
            obj.L = pinv(G * F);
            
            % errors will have the initial point too
            obj.L = [zeros(N*dim_u,dim_y),obj.L];
        end
        
        function formL(obj,model,traj)
            
            dim = model.SIM.dimx/2; % dof for positions
            K = traj.K;
            N = traj.N - 1;
            
            % construct learning matrix first
            %obj.basicL(traj);
            obj.lift(model,traj);
            
            if obj.modifyPosOnly
                Kl = zeros(dim*N);
                % construct lifted domain matrix Kl
                for l = 1:N
                    vec = (l-1)*dim + 1:l*dim;
                    Kl(vec,vec) = K(1:dim,1:dim,l); % on diagonals only
                end
            else
                % take pinv(K) to reflect back on the trajectories s            
                Kl = zeros(size(K,1)*N,size(K,2)*N);
                % construct lifted domain matrix Kl
                for l = 1:N
                    vec1 = (l-1)*size(K,1) + 1:l*size(K,1);
                    vec2 = (l-1)*size(K,2) + 1:l*size(K,2);
                    Kl(vec1,vec2) = K(:,:,l); % on diagonals only
                end
            end

            obj.L = pinv(Kl)*obj.L;
            %obj.L = Kl \ obj.L;
                
            
        end
        
        %% basis functions are unscaled gaussians
        function out = basis(obj,t,h,c)
            out = exp(-h * (t - c).^2);
        end
        
        %% update weights of the trajectory or the trajectory directly
        function traj2 = updateWeights(obj,traj,y)
            
            h = traj.t(2) - traj.t(1);
            N = traj.N - 1;
            K = traj.K;
            Cout = obj.C;
            dim = size(Cout,2);
            dev = y - traj.s;   
            
            % ILC update
            w = obj.inp_last + obj.Psi \ obj.L*dev(:);
            obj.inp_last = w;
            % form back the trajectory
            sbar = obj.Psi * w;
            sbar = reshape(sbar,dim,N);
            %s = [Cout*sbar,traj.s(:,end)];
            s = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t, s, traj.unom, K);
            
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
            Sv = repmat(eye(dim),N,N);
            S = tril(Sv);
            
            % ILC update
            wd = obj.inp_last + pinv(obj.Psi)*D*obj.L*dev(:);
            obj.inp_last = wd;
            % form back the trajectory
            sbar = S * obj.Psi * wd;
            sbar = reshape(sbar,dim,N);
            %s = [Cout*sbar,traj.s(:,end)];
            s = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t, s, traj.unom, K);
            
            % get s from derivatives
            %S = h*tril(ones(N+1));
            %sbar = S * sd;           

            
        end
        
        % update trajectories 
        function traj2 = updateTraj(obj,traj,y)
            
            N = traj.N - 1;
            K = traj.K;
            dev = y - traj.s;
            %dev = dev(:,2:end);
            Cout = obj.C;
            dim = size(Cout,2);

            % ILC update
            if obj.modifyPosOnly
                posdiff = obj.L * dev(:);
                % add zeros
                posdiff = reshape(posdiff,dim/2,N);
                state = reshape(obj.inp_last,dim,N);
                sbar = state + [posdiff;zeros(dim/2,N)];
                obj.inp_last = sbar(:);
            else
                s = obj.inp_last + obj.L*dev(:);
                obj.inp_last = s;
                % form back the trajectory
                sbar = reshape(s,dim,N);
            end
            %s = [Cout*sbar,traj.s(:,end)];
            s = [sbar,sbar(:,end)];
            
            traj2 = Trajectory(traj.t,s,traj.unom,K);
            
        end
        
        %% update DMPs
        function dmp = updateDMP(obj,dmp,traj,y)
            
            N = traj.N - 1;
            h = traj.t(2) - traj.t(1);
            ref = traj.s;
            dev = y - ref;
            Cout = obj.C;
            dim = size(Cout,2);  
            
            s = obj.inp_last + obj.L*dev(:);
            obj.inp_last = s;
            % form back the trajectory
            sbar = reshape(s,dim,N);
            s = Cout*[sbar,sbar(:,end)];            
            
            q = s;
            qd = diff(q')'/h;
            qd(:,end+1) = qd(:,end);
            qdd = diff(qd')'./h;
            qdd(:,end+1) = qdd(:,end);
    
            goals = q(:,end);
            y0 = [q(:,1),qd(:,1)];

            for i = 1:length(dmp)

                % set/update goal and initial state
                dmp(i).setInitState(y0(i,:)');
                % initial states of DMPs
                dmp(i).setGoal(goals(i));
                %dmp(i) = discreteDMP(can,alpha,beta,goal(1),yin,numbf);
                dmp(i).regressLive(q(i,:)',qd(i,:)',qdd(i,:)',goals(i));

            end
            
            %Fs = dmp.FORCE.Fs;
            %w_next = obj.inp_last - pinv(Fs)*obj.L*dev(:);
            %dmp.FORCE.w = w_next;
            %obj.inp_last = w_next;
            
        end
        
        %% switching function that switches between different
        % methods. In the first three, out is a new trajectory class
        % updateDMP outputs a DMP class
        function out = feedforward(obj,traj,dmps,y)
        
            update_method = obj.upd_meth;
            
            switch update_method 
                case 't'
                    out = obj.updateTraj(traj,y);
                case 'w'
                    out = obj.updateWeights(traj,y);
                case 'wd'
                    out = obj.updateDerivativeWeights(traj,y);
                otherwise %updates dmp weights
                    out = obj.updateDMP(dmps,traj,y);
                    
            end
            
        end
    end
    
end