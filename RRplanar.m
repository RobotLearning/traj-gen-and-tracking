% Two-link planar revolute robot manipulator (RR)

classdef RRplanar < Robot

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % cost function structure (handle and weight matrix)
        COST
        % input bound
        bound
        % cartesian coordinates of the trajectory
        x, xd, xdd
        % joint space coordinates
        q, qd, qdd
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            % initialize everything to zero
            obj.PAR.const.g = 0;
            obj.PAR.link1.mass = 0;
            obj.PAR.link2.mass = 0;
            obj.PAR.link1.length = 0;
            obj.PAR.link2.length = 0;
            obj.PAR.link1.centre.dist = 0;
            obj.PAR.link2.centre.dist = 0;
            obj.PAR.link1.inertia = 0;
            obj.PAR.link2.inertia = 0;
            obj.PAR.link1.motor.inertia = 0;
            obj.PAR.link2.motor.inertia = 0;
            obj.PAR.link1.motor.gear_ratio = 0;
            obj.PAR.link2.motor.gear_ratio = 0;
                         
            % check that the input has all the fields
            % TODO: is there a better way?
            assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            % initialize all fields
            obj.CON.link1.q.max = Inf;
            obj.CON.link1.q.min = -Inf;
            obj.CON.link1.u.max = Inf;
            obj.CON.link1.u.min = -Inf;
            obj.CON.link2.q.max = Inf;
            obj.CON.link2.q.min = -Inf;
            obj.CON.link2.u.max = Inf;
            obj.CON.link2.u.min = -Inf;
            
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % copies the rectangular bounds on R2 input space
        function set.bound(obj, mat)
            
            if isempty(obj.bound)
                obj.bound = Inf * [-1,1;
                                   -1,1];
            end
            
            % make sure the values are tighter
            obj.bound(:,1) = max(obj.bound(:,1), mat(:,1));
            obj.bound(:,2) = min(obj.bound(:,2), mat(:,2));
            
        end
        
        % change the cost function
        function set.COST(obj, Q)
            obj.COST.Q = Q;
            obj.COST.fnc = @(x1,x2) (x1-x2)'*Q*(x1-x2);
        end
        
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = RRplanar(par,con,cost)
            
            % simulation time step
            obj.SIM.h = 0.02;
            % dimension of the x vector
            obj.SIM.dim_x = 2;
            % dimension of the action space
            obj.SIM.dim_u = 2;
            % noise variance of the model (white noise)
            obj.SIM.eps = 0.00003;
            
            % set object parameter
            obj.PAR = par;
            % set object disturbance
            obj.DIST = dist;
            % set object constraints
            obj.CON = con;
            
            % bounds of the inputs
            bnd = [con.fmin, con.fmax; 
                  -con.phi_dot_max, con.phi_dot_max]; 
            % set object bounds
            obj.bound = bnd;
                     
            % cost function handle
            obj.COST = cost.Q;
        end
        
        % provides nominal model
        function [x_dot,varargout] = nominal(~,obj,x,u,flg)
            % differential equation of the inverse dynamics
            % x_dot = Ax + B(x)u + C
            x_dot = inverseDynamics(x,u);
            
            % return jacobian matrices
            if flg
                vec = [0; -u(1)*cos(x(5)); 0; u(1)*sin(x(5)); 0];
                dfdx = [A(:,1:4), vec];
                varargout{1} = dfdx;
                varargout{2} = B;
            end
        end
        
        % provides actual model
        function x_dot = actual(~,obj,x,u,flg)
            
            % TODO: change parameters
            
            % differential equation of the inverse dynamics
            % x_dot = Ax + B(x)u + C
            x_dot = inverseDynamics(x,u);
            
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift(obj,trj)
            
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.dim_x;
            dimu = obj.dim_u;
            umin(1,:) = obj.CON.fmin - u_trj(1,:);
            umin(2,:) = -obj.CON.phi_dot_max - u_trj(2,:);
            umax(1,:) = obj.CON.fmax - u_trj(1,:);
            umax(2,:) = obj.CON.phi_dot_max - u_trj(2,:);

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

        % wrapper for the splines-based trajectory generator method
        % using the differential flatness of the dynamics
        function Traj = trajectory(obj,shape,spline)

            spline_x = spline(1,:);
            spline_y = spline(2,:);
            [t,u_nom,x_nom] = quad_traj_gen(shape,obj.PAR,obj.CON,...
                                            obj.h,spline_x,spline_y);
            fun0 = @(t,obj,x,u) nominal(t,obj,x,u,false);
            fun = @(t,obj,x,u) nominal(t,obj,x,u,false) ...
                             + disturbance(t,obj,x,u);
            x_real = simulate(obj,t,x_nom(:,1),u_nom,fun);
            x_pred = predict_full(obj,t,x_real,u_nom,fun0);
            Traj = Trajectory(t,spline,x_nom,u_nom);
            cost = obj.COST;
            % trajectory generator generates an extra input
            Traj.addPerformance(u_nom(:,1:end-1),x_real,...
                                x_pred,cost,'splines');
            % TODO add noise to x! (as observed variable)
        end
        
    end
end