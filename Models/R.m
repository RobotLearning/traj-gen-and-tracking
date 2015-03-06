% One-link revolute arm (pendulum)

classdef R < Robot

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
        % jacobian matrix
        jac
        % learning in joint space?
        flag_jspace
        % represent references in joint space
        flag_ref_jsp
    end
    
    methods
        
        %% Constructor and initializing methods 
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            % initialize everything to zero
            obj.PAR.const.g = 0;
            obj.PAR.link.mass = 0;
            obj.PAR.link.length = 0;
            obj.PAR.link.centre.dist = 0;
            obj.PAR.link.inertia = 0;
            obj.PAR.link.motor.inertia = 0;
            obj.PAR.link.motor.gear_ratio = 0;
                         
            % check that the input has all the fields
            % TODO: is there a better way?
            %assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
            
            % set observation matrix
            obj.C = STR.C;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            % initialize all fields
            obj.CON.link.q.max = Inf;
            obj.CON.link.q.min = -Inf;
            obj.CON.link.qd.max = Inf;
            obj.CON.link.qd.min = -Inf;
            obj.CON.link.qdd.max = Inf;
            obj.CON.link.qdd.min = -Inf;
            obj.CON.link.u.max = Inf;
            obj.CON.link.u.min = -Inf;
            obj.CON.link.udot.max = Inf;
            obj.CON.link.udot.min = -Inf;
            
            % check that the input has all the fields
            %assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
            obj.SIM.eps_d = sim.eps_d;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
            % TODO: matlab is right to complain here
            obj.flag_jspace = ~sim.cartesian;
            obj.flag_ref_jsp = sim.jref;
        end
        
        % change the cost function
        function set.COST(obj, cost)
            obj.COST.Q = cost.Q;
            obj.COST.R = cost.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*cost.Q*(x1-x2));
            %assert(length(cost.Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = R(par,con,cost,sim)
            
            obj.SIM = sim;            
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;        
            % cost function handle
            obj.COST = cost;
            % TODO: construct jacobian
            obj.jac = [];
        end
        
        %% Methods for evolving dynamics models
        
        % provides nominal model
        function [x_dot,varargout] = nominal(obj,~,x,u,flg)
            % differential equation of the inverse dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            [x_dot,dfdx,dfdu] = obj.dynamics(x,u,true);
            
            % return jacobian matrices
            if flg
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            end
        end
        
        % provides actual model
        % TODO: should we wrap the dynamics?
        function x_dot = actual(obj,~,x,u)
            
            % change parameters
            par.const.g = obj.PAR.const.g;
            par.link.mass = 1.0 * obj.PAR.link.mass;
            par.link.length = 1.0 * obj.PAR.link.length;
            par.link.centre.dist = obj.PAR.link.centre.dist;
            par.link.inertia = 1.0 * obj.PAR.link.inertia;            
            par.link.motor.inertia = 1.0 * obj.PAR.link.motor.inertia;            
            par.link.motor.gear_ratio = 1.0 * obj.PAR.link.motor.gear_ratio;            
            
            % differential equation of the inverse dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            x_dot = RDynamics(x,u,par,false);
            
        end
        
        %% Kinematics and dynamics models
        
        % run kinematics
        function x = kinematics(obj,q)
            
            l1 = obj.PAR.link.length;
            x = [l1 * cos(q); l1 * sin(q)];
        end
        
        % call inverse kinematics from outside
        function q = invKinematics(obj,x)
            
            q = atan(x(2,:)./x(1,:));
        end
                   
        % dynamics to get u
        function u = invDynamics(obj,q,qd,qdd)
            u = RInvDynamics(q,qd,qdd,obj.PAR);
        end
        
        % dynamics to qet Qd = [qd,qdd]
        function [Qd, varargout] = dynamics(obj,Q,u,flag)
            if flag
                [Qd, dfdx, dfdu] = RDynamics(Q,u,obj.PAR,flag);
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            else
                Qd = RDynamics(Q,u,obj.PAR,flag);
            end
        end
        
        %% Make an animation of the robot manipulator
        function animateArm(obj,q_actual,s)
            x = obj.kinematics(q_actual);
            animateR(x,s);
        end
        
        % unused except for aILC
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            h = obj.SIM.h;
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            umin = obj.CON.link.u.min - u_trj;
            umax = obj.CON.link.u.max - u_trj;

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/h;
            D = D(1:end-dimu,:);
            
            u_dot_max = obj.CON.link.udot.max;
            u_dot_min = obj.CON.link1.udot.min;
            U_dot_max = repmat(u_dot_max,N-1,1);
            U_dot_min = repmat(u_dot_min,N-1,1);
            u_star = u_trj(:);

            L = [D; -D];
            q = [U_dot_max - D*u_star;
                 -U_dot_min + D*u_star];
            
        end
        
    end
end