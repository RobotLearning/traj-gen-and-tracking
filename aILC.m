% Iterative Learning Control using Angela's Kalman-filtering
% approach

classdef aILC < ILC
    
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
        % L and q are (convex) objective function penalties
        L
        q
        % scaling matrix 
        S
        % scales for the covariances (process/output)
        eps
        eps_d
        % Kalman filter
        filter
        
    end
    
    methods
        
        function obj = aILC(model,trj)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'Kalman filter based ILC';
            obj.error = 0;
            
            N = trj.N - 1;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;

            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = eye(N*dim_x);
            obj.H = zeros(N*dim_x,N*dim_u); 
            obj.eps = model.SIM.eps;
            obj.eps_d = model.SIM.eps_d;
            obj.u_last = zeros(dim_u,N);
            
            obj.lift(model,trj);
            % initialize Kalman filter
            mats.M = obj.eps * eye(N*dim_x); % covariance of process x measurement noise
            mats.Omega = obj.eps_d * eye(N*dim_x); % initial covariance of d noise
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
                    vec_x = (l-1)*dim_x + 1:l*dim_x;
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
            [obj.umin,obj.umax,obj.L,obj.q] = model.lift_constraints(trj);
            
            % construct scaling matrix S
            Sx = eye(dim_x);
            % weighing trajectory deviation by these values
            Sw = sqrt(model.COST.Q);
            Tw = eye(dim_x);
            obj.S = cell(1,N);
            [obj.S{:}] = deal(Tw * Sx * Sw);
            obj.S = blkdiag(obj.S{:});

        end
        
        function us = feedforward(obj,trj,model,dev)
            
            h = model.SIM.h;
            % get rid of x0 in dev
            dev = dev(:,2:end);            
            N = trj.N - 1;
            %dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            F = obj.F;
            S = obj.S;
            L = obj.L;
            q = obj.q;
            umin = obj.umin; 
            umax = obj.umax;
            
            % input deviation penalty matrix D
            D0 = eye(dim_u*N); 
            D1 = (diag(ones(1,dim_u*(N-1)),dim_u) - eye(dim_u*N))/h;
            D1 = D1(1:end-dim_u,:); % D1 is (N-1)*dimu x N*dimu dimensional
            D2 = (diag(ones(1,dim_u*(N-2)),2*dim_u) - ... 
                 2*diag(ones(1,dim_u*(N-1)),dim_u) + eye(dim_u*N)) ...
                 /(h^2);
            D2 = D2(1:end-2*dim_u,:); % D2 is (N-2)*dimu x N*dimu dimensional
            % penalty scale
            a0 = 5e-5;
            a1 = 5e-5; 
            a2 = 5e-5;
            % slack for input u
            %eps_u = 1e-6;
            
            a = a2;
            D = D2;
            
            % get disturbance estimate from filter
            u_pre = obj.u_last(:);
            obj.filter.predict(u_pre);
            obj.filter.update(dev,u_pre);
            d = obj.filter.x;
    
            % solve with quadprog
            options = optimset('Display', 'iter', ...
                      'Algorithm', 'interior-point-convex');
            uiter = quadprog(2*(F'*S')*S*F + 2*a*(D'*D), ...
                          2*F'*S'*d, L, q, [], [], umin, umax, [], options);
            
            obj.u_last = reshape(uiter,dim_u,N);
            us = trj.unom(:,1:N) + obj.u_last;
        end
        
    end
    
end