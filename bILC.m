% Basic Iterative Learning Control using PD-type input update

classdef bILC < ILC
    
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
        
        % ILC's Last input sequence
        u_last
        % Lifted state matrix F
        F
        % Lifted state matrix G
        G
        % Lifted state matrix H
        H
        % Lifted state inequalities
        umin
        umax
    end
    
    % fields that are particular to this implementation
    properties
        
        % scales for the covariances (process/output)
        eps
        eps_M
        % Kalman filter
        filter
        
    end
    
    methods
        
        function obj = bILC(model,trj)
                        
            obj.episode = 0;
            obj.color = 'b';
            obj.name = 'Arimoto type ILC';
            obj.error = 0;
            
            N = trj.N - 1;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;

            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = eye(N*dim_x);
            obj.H = zeros(N*dim_x,N*dim_u); 
            obj.eps = 0.003;
            obj.eps_M = 0.05;
            obj.u_last = trj.unom(:,1:N);
            
            obj.lift(model,trj);
            % initialize Kalman filter
            mats.M = obj.eps_M * eye(N*dim_x); % covariance of process x measurement noise
            mats.Omega = obj.eps * eye(N*dim_x); % initial covariance of d noise
            mats.A = eye(N*dim_x);
            mats.B = zeros(N*dim_x,N*dim_u);
            mats.C = obj.G;
            mats.D = obj.G * obj.F + obj.H;
            obj.filter = Filter(model,trj,mats);
            
        end
        
        % get the lifted vector representation 
        % around the trajectory
        % TODO: better way to construct F?
        function obj = lift(obj,model,trj)
            
            N = trj.N - 1;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            
            % get linear time variant matrices around trajectory
            [A,B] = model.linearize(trj);
            
            % construct lifted domain matrix F
            for l = 1:N
                for m = 1:N
                    vec_x = (l-1)*dim_x + 1: l*dim_x;
                    vec_u = (m-1)*dim_u + 1:m*dim_u;
                    if m <= l
                        mat = eye(dim_x);
                        for i = m+1:l
                            mat = mat * A(:,:,i);
                        end
                        obj.F(vec_x,vec_u) = mat * B(:,:,m); % on diagonals only B(:,m)
                    end
                end
            end
            
            % construct umin and umax
            % extract constraints
            %[obj.umin,obj.umax,obj.L,obj.q] = model.lift_constraints(trj);
            
            % construct scaling matrix S
            %Sx = eye(dim_x);
            % weighing trajectory deviation by these values
            %Sw = sqrt(model.COST.Q);
            %Tw = eye(dim_x);
            %obj.S = cell(1,N);
            %[obj.S{:}] = deal(Tw * Sx * Sw);
            %obj.S = blkdiag(obj.S{:});

        end
        
        function u_next = feedforward(obj,trj,model,dev)
            
            %h = model.SIM.h;
            % get rid of x0 in dev
            ddev = diff(dev(1:2,:)')';
            dev = dev(1:2,2:end);            
            %N = trj.N - 1;
            %dim_x = model.SIM.dimx;
            %dim_u = model.SIM.dimu;
            %F = obj.F;
            %S = obj.S;
            %L = obj.L;
            %q = obj.q;
            %umin = obj.umin; 
            %umax = obj.umax;
            
            % get disturbance estimate from filter
            %u_pre = obj.u_last(:);
            %obj.filter.predict(u_pre);
            %obj.filter.update(dev,u_pre);
            %d = obj.filter.x;
    
            % set learning rate
            alpha_p = 0.0 * max(max(obj.u_last))/max(max(abs(dev)));
            alpha_d = 0.2 * max(max(obj.u_last))/max(max(abs(ddev)));
            u_next = obj.u_last + alpha_p * dev + alpha_d * ddev;
            
            obj.u_last = u_next;
            
        end
        
    end
    
end