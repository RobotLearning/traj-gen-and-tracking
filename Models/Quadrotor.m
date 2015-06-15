% Quadrotor 2D-model implementing the Model superclass

classdef Quadrotor < Model

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % cost function structure (handle and weight matrix)
        COST
        % fields necessary for simulation and plotting, noise etc.
        SIM
        % observation matrix
        C        
    end
    
    methods
        
        %% Methods for initializing variables and constructor
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            % initalize all fields to 0
            obj.PAR = struct('Iy',0,'L',0,'m',0,'g',0);
            obj.PAR.Quad = struct('A',0);
                         
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            % initialize all fields to 0
            obj.CON = struct('fi_max',0,'fi_min',0,'fi_dot_max',0,...
                             'fmax',0,'fmin',0,'f_dot_max',0, 'phi_max',0,...
                             'phi_dot_max',0,'phi_ddot_max',0,...
                             'NumT',0);
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.CON), fieldnames(obj.CON))));
            obj.CON = STR;
            
        end
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            obj.SIM.eps_m = sim.eps_m;
            %assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
            %       'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        function set.COST(obj, cost)
            obj.COST.Q = cost.Q;
            obj.COST.R = cost.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*cost.Q*(x1-x2));
            %assert(length(Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        function obj = Quadrotor(par,con,cost,sim)
            
            obj.SIM = sim;            
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;        
            % cost function handle
            obj.COST = cost;
            
            % bounds of the inputs
            %bnd = [con.fmin, con.fmax; 
            %      -con.phi_dot_max, con.phi_dot_max]; 
            % set object bounds
            %obj.bound = bnd;
        end
        
        %% Methods for evolving dynamics models
        % provides nominal model
        function [x_dot,varargout] = nominal(obj,t,x,u,flg)
            % x_dot = quadrocopterNominalDynamics(t,x,u)
            % differential equation of the quadrocopter
            % Dynamics
            % x_dot = Ax + B(x)u + C
            A = [0 1 0 0 0;
                 zeros(1,5);
                 0 0 0 1 0;
                 zeros(2,5)];
            B = [0          0;
                 -sin(x(5)) 0;
                 0          0;
                 cos(x(5))  0;
                 0          1];
            C = [0; 0; 0; -obj.PAR.g; 0];
            x_dot = A*x + B*u + C;
            
            % return jacobian matrices
            if flg
                vec = [0; -u(1)*cos(x(5)); 0; u(1)*sin(x(5)); 0];
                dfdx = [A(:,1:4), vec];
                varargout{1} = dfdx;
                varargout{2} = B;
            end
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % add disturbance
            x_dot = obj.nominal(t,x,u,false) + obj.disturbance(t,x,u);            
            
        end
        
        % provides disturbance to be added to the nominal model
        % TODO: DIST structure needs to be replaced
        function delta = disturbance(obj,t,x,u)
             
            type = 'gravity';
            
            switch type
                case 'gravity'
                    F_z = obj.PAR.g - obj.DIST.g;
                    delta = [0; 0; 0; F_z; 0];                
                case 'wind'
                    F_y = obj.DIST.BlowPressure * obj.PAR.Quad.A * ...
                          sin(obj.DIST.Angle + x(5)) * cos(obj.DIST.Angle);
                    F_z = obj.DIST.BlowPressure * obj.PAR.Quad.A * ...
                          sin(obj.DIST.Angle + x(5)) * sin(obj.DIST.Angle);
                    delta = [0; F_y; 0; F_z; 0];
                case 'actuator'
                    if ~isfield(obj.DIST, 'k')
                        obj.DIST.k = 0.1; % multiply B matrix by 1+k
                    end
                    delta_B = [0 0;
                              -obj.DIST.k*sin(x(5)) 0;
                               0 0;
                               obj.DIST.k*cos(x(5)) 0;
                               0 1];
                    delta = delta_B * u;
            end
        end        
        
        %% get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
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

        %% splines-based trajectory generator method
        % taking advantage of the differential flatness of the dynamics
        function Traj = generateInputs(obj,shape,spline)

            spline_x = spline(1,:);
            spline_y = spline(2,:);
            [t,u_nom,x_nom] = quad_traj_gen(shape,obj.PAR,obj.CON,...
                                            obj.h,spline_x,spline_y);
            Traj = Trajectory(t,spline,u_nom,[]);
        end
        
    end
end