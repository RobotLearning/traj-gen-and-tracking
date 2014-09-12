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
        eps_M
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
            dim_x = model.dim_x;
            dim_u = model.dim_u;
            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = eye(N*dim_x);
            obj.H = zeros(N*dim_x,N*dim_u); 
            obj.eps = 0.003;
            obj.eps_M = 0.05;
            obj.u_pre = zeros(dim_u*N,1);
            
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
            dim_x = model.dim_x;
            dim_u = model.dim_u;
            
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
            [obj.umin,obj.umax,obj.L,obj.q] = model.lift(trj);
            
            % construct scaling matrix S
            Sx = eye(dim_x);
            % weighing trajectory deviation by these values
            Sw = sqrt(model.COST.Q);
            Tw = eye(dim_x);
            obj.S = cell(1,N);
            [obj.S{:}] = deal(Tw * Sx * Sw);
            obj.S = blkdiag(obj.S{:});

        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj)
            
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.dim_x;
            dimu = obj.dim_u;
            umin(1,:) = obj.CON.link1.u.min - u_trj(1,:);
            umin(2,:) = obj.CON.link2.u.min - u_trj(2,:);
            umax(1,:) = obj.CON.link1.u.max - u_trj(1,:);
            umax(2,:) = obj.CON.link2.u.max - u_trj(2,:);

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/obj.h;
            D = D(1:end-dimu,:);
            % construct L1 and L2
            L1 = zeros(N-1, N*dimu);
            L2 = zeros(N-1, N*dimu);
            a = 1/4;
            b = obj.PAR.Iy/(2 * obj.PAR.m * obj.PAR.L); 
            b = b/obj.h; %b_bar
            vec1 = [a -b 0 b];
            vec2 = [a b 0 -b];
            for i = 1:N-1
                L1(i,:) = [zeros(1,(i-1)*dimu), vec1, zeros(1,(N-i-1)*dimu)];
                L2(i,:) = [zeros(1,(i-1)*dimu), vec2, zeros(1,(N-i-1)*dimu)];
            end
            u_dot_max = [4*obj.CON.fi_dot_max; obj.CON.phi_ddot_max];
            U_dot_max = repmat(u_dot_max,N-1,1);
            u_star = u_trj(:); 

            L = [D; -D; L1; -L1; L2; -L2];
            q = [U_dot_max - D*u_star; 
                 U_dot_max + D*u_star;
                 obj.CON.fmax - L1*u_star;
                 -obj.CON.fmin + L1*u_star;
                 obj.CON.fmax - L2*u_star;
                 -obj.CON.fmin + L2*u_star];    
            
        end
        
        function us = feedforward(obj,trj,model,dev)
            
            % get rid of x0 in dev
            dev = dev(:,2:end);            
            N = trj.N - 1;
            %dim_x = model.dim_x;
            dim_u = model.dim_u;
            
            % input deviation penalty matrix D
            %D0 = eye(length(u_last)); 
            %D1 = (diag(ones(1,dim_u*(N-1)),dim_u) - eye(dim_u*N))/model.h;
            %D1 = D1(1:end-dim,:); % D1 is (N-1)*dimu x N*dimu dimensional
            D2 = (diag(ones(1,dim_u*(N-2)),2*dim_u) - ... 
                 2*diag(ones(1,dim_u*(N-1)),dim_u) + eye(dim_u*N)) ...
                 /(model.h^2);
            D2 = D2(1:end-2*dim_u,:); % D2 is (N-2)*dimu x N*dimu dimensional
            % penalty scale
            %a0 = 5e-5; 
            %a1 = 5e-5; 
            a2 = 5e-5;
            % slack for input u
            %eps_u = 1e-6;
            
            % get disturbance estimate from filter
            obj.filter.predict(obj.u_last);
            obj.filter.update(dev,obj.u_last);
            d = obj.filter.x;
    
            % solve with quadprog
            options = optimset('Display', 'iter', ...
                      'Algorithm', 'interior-point-convex');
            uiter = quadprog(2*(obj.F'*obj.S')*obj.S*obj.F + 2*a2*(D2'*D2), ...
                          2*obj.F'*obj.S'*d, obj.L, obj.q, [], [], ...
                          obj.umin, obj.umax, [], options);
            
            us = trj.unom(:,1:N) + reshape(uiter,dim_u,N);
        end
        
    end
    
end