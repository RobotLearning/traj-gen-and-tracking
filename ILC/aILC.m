% Iterative Learning Control using Kalman-filtering
% approach (see Angela Schoellig's paper)
%

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
        inp_last
        % Lifted state matrix F
        F
        % Lifted state matrix G
        G
        % Lifted state matrix H
        H
    end
    
    % fields that are particular to this implementation
    properties
        % scaling matrix 
        S
        % scale for the process/output covariance
        eps
        % scale for the disturbance (d) covariance 
        eps_d
        % Kalman filter
        filter
        
    end
    
    methods
        
        %% Constructor for convex optim + Kalman filter based ILC
        function obj = aILC(model,trj)
                        
            obj.episode = 0;
            obj.color = 'k';
            obj.name = 'Kalman filter based ILC';
            obj.error = 0;
            
            N = trj.N - 1;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            dim_y = model.SIM.dimy;

            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = eye(N*dim_y, N*dim_x);
            obj.H = zeros(N*dim_y,N*dim_u); 
            obj.eps = model.SIM.eps;
            obj.eps_d = max(1e-6,model.SIM.eps);
            obj.inp_last = trj.unom(:,1:N);
            
            obj.lift(model,trj);
            % initialize Kalman filter
            mats.M = obj.eps * eye(N*dim_y); % covariance of process x measurement noise
            mats.O = obj.eps_d * eye(N*dim_x); % initial covariance of d noise
            mats.A = eye(N*dim_x);
            mats.B = zeros(N*dim_x,N*dim_u);
            mats.C = obj.G;
            mats.D = obj.G * obj.F + obj.H;
            obj.filter = Filter(model,trj,mats);
            
        end
        
        %% get the lifted vector representation 
        % around the trajectory
        function obj = lift(obj,model,trj)
            
            N = trj.N - 1;
            dim_y = model.SIM.dimy;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            
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
                        obj.F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                    end
                end
                
            else
            
                % get linear time variant matrices around trajectory
                [Ad,Bd] = model.linearize(trj);

                % construct lifted domain matrix F
                for l = 1:N
                    for m = 1:N
                        vec_x = (l-1)*dim_x + 1:l*dim_x;
                        vec_u = (m-1)*dim_u + 1:m*dim_u;
                        if m <= l
                            mat = eye(dim_x);
                            for i = m+1:l
                                mat = mat * Ad(:,:,i);
                            end
                            obj.F(vec_x,vec_u) = mat * Bd(:,:,m); % on diagonals only B(:,m)
                        end
                    end
                end
            end
            
            % construct scaling matrix S
            Sx = eye(dim_x);
            % weighing trajectory deviation by these values
            % all nonoutput states should go to zero
            Q2 = diag([zeros(1,dim_y), ones(1,dim_x-dim_y)]);
            Cbar = eye(dim_x) - model.C'*((model.C*model.C')\model.C);
            Sw = Cbar'*Q2*Cbar + model.C'*model.COST.Q*model.C;
            Tw = eye(dim_x);
            obj.S = cell(1,N);
            % Sx: scaling matrix that scales different state dimensions
            % Sw: weight different state dimensions by importance
            % Tw: change the scaling along the trajectory
            [obj.S{:}] = deal(Tw * Sx * Sw);
            obj.S = blkdiag(obj.S{:});

        end
        
        %% Main ILC function
        function unext = feedforward(obj,traj,model,y)
                        
            ydev = y - traj.s;
            % deviation vector starts from x(1)
            dev = ydev(:,2:end);
            
            Nu = traj.N - 1;            
            dim_u = model.SIM.dimu;
            h = model.SIM.h;
            % get u from last applied input
            u = obj.inp_last - traj.unom(:,1:Nu);
            
            % Filter to estimate disturbance
            obj.filter.predict(u(:));
            obj.filter.update(dev(:),u(:));
            d = obj.filter.x;

            % Constraints and penalty matrices
            % construct umin and umax
            % extract constraints
            % arrange in format Lu <= q
            % call the particular constraints-generating code
            if ismethod(model,'lift_constraints')
                [umin,umax,L,q] = model.lift_constraints(traj,obj);
            else
                umin = []; umax = [];
            end           

    
            % input deviation penalty matrix D
            D0 = eye(dim_u*Nu); 
            D1 = (diag(ones(1,dim_u*(Nu-1)),dim_u) - eye(dim_u*Nu))/h;
            D1 = D1(1:end-dim_u,:); % D1 is (Nu-1)*nu x Nu*nu dimensional
            D2 = (diag(ones(1,dim_u*(Nu-2)),2*dim_u) - ... 
                 2*diag(ones(1,dim_u*(Nu-1)),dim_u) + eye(dim_u*Nu))/(h^2);
            D2 = D2(1:end-2*dim_u,:); % D2 is (Nu-2)*nu x Nu*nu dimensional
            % penalty scale
            a0 = 5e-5; a1 = 5e-5; a2 = 5e-5;
    
            % solve with quadprog
            options = optimset('Display', 'iter', 'Algorithm', 'interior-point-convex');
            M0 = obj.S*obj.F;
            M = (M0'*M0) + a2*(D2'*D2);
            v = obj.F'*obj.S'*d;
            u = quadprog(M, v, [], [], [], [], umin, umax, [], options);
            
            unext = traj.unom(:,1:Nu) + reshape(u,dim_u,Nu);
           
            
        end
        
    end
    
end