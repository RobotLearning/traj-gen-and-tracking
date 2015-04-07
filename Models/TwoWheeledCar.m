% Two-wheeled 2D car or robot kinematical model 

classdef TwoWheeledCar < Model

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
            
            % initialize everything to zero
            obj.PAR.wheel1.radius = 0;
            obj.PAR.wheel2.radius = 0;
            obj.PAR.length = 0;
                         
            % TODO: check that the input has all the fields
            obj.PAR = STR;
            
            % set observation matrix
            obj.C = STR.C;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % initialize all fields to zero
            obj.CON.state.x.max = 0;
            obj.CON.state.x.min = 0;
            obj.CON.state.y.max = 0;
            obj.CON.state.y.min = 0;
            obj.CON.wheel1.u.max = 0;
            obj.CON.wheel1.u.min = 0;
            obj.CON.wheel2.u.max = 0;
            obj.CON.wheel2.u.min = 0;
            obj.CON.wheel1.udot.max = 0;
            obj.CON.wheel1.udot.min = 0;
            obj.CON.wheel2.udot.max = 0;
            obj.CON.wheel2.udot.min = 0;
            
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        function set.COST(obj, STR)
            obj.COST.Q = STR.Q;
            obj.COST.R = STR.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*STR.Q*(x1-x2));
            %assert(length(Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = TwoWheeledCar(par,con,cost,sim)
            
            obj.SIM = sim;
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;    
            % cost function handle
            obj.COST = cost;
        end
        
        %% Methods for evolving dynamics models
        
        % provides nominal model
        function [x_dot,varargout] = nominal(obj,t,x,u,flg)
            % differential equation of the inverse dynamics
            % x_dot = B(x)u
            [x_dot,dfdx,dfdu] = obj.differentialKinematics(t,x,u,true);
            
            % return jacobian matrices
            if flg
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            end
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % change parameters
            par.wheel1.radius = 1.1 * obj.PAR.wheel1.radius;
            par.wheel2.radius = 0.9 * obj.PAR.wheel2.radius;
            par.length = 1.0 * obj.PAR.length;
            % time varying drift on the angle
            %a = 1*[0;0;1];
            
            % differential equation of the inverse dynamics
            x_dot = robotTwoWheelsKinematics(t,x,u,par,false);
            %x_dot = robotTwoWheelsKinematics(t,x,u,par,false) + a*sin(2*pi*t);
            
        end
        
        %% Kinematics based simple model
                           
        % dynamics to get xdot
        function [x_dot,varargout] = differentialKinematics(obj,t,x,u,flag)
            if flag
                [x_dot,dfdx,dfdu] = robotTwoWheelsKinematics(t,x,u,obj.PAR,true);
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            else
                x_dot = robotTwoWheelsKinematics(t,x,u,obj.PAR,false);
            end
        end
        
        %% Make an animation of the model
        function animate(obj,x,s)
            
            animateCar(x,s,obj.PAR);
        end

        %% Generate inputs for tracking reference
        % using a simple inverse kinematics method
        % TODO: extend using planning to incorporate constraints
        function Traj = generateInputs(obj,t,x_des)

            N = length(t)-1;
            dimu = obj.SIM.dimu;
            h = obj.SIM.h;
            unom = zeros(dimu,N);
            R1 = obj.PAR.wheel1.radius;
            R2 = obj.PAR.wheel2.radius;
            d = obj.PAR.length;
            Cout = obj.C;

            for i = 1:N
                % discretize B(phi) around phi0
                B = h * [R1/2 * cos(x_des(3,i)), R2/2 * cos(x_des(3,i));
                         R1/2 * sin(x_des(3,i)), R2/2 * sin(x_des(3,i));
                         R1/d,                   -R2/d];
                % deviation along desired trajectory
                delta = x_des(2:3,i+1) - x_des(2:3,i);
                % solve for u
                unom(:,i) = B(2:3,:)\delta;
            end
            
            Traj = Trajectory(t,Cout*x_des,unom,[]);
        end
        
        % unused except for aILC
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            h = obj.SIM.h;
            N = trj.N - 1; 
            s = trj.s;
            u_trj = trj.unom(:,1:N);
            %dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            
            % input constraints
            umin(1,:) = obj.CON.wheel1.u.min - u_trj(1,:);
            umin(2,:) = obj.CON.wheel2.u.min - u_trj(2,:);
            umax(1,:) = obj.CON.wheel1.u.max - u_trj(1,:);
            umax(2,:) = obj.CON.wheel2.u.max - u_trj(2,:);

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/h;
            D = D(1:end-dimu,:);
            
            u_dot_max = [obj.CON.wheel1.udot.max; obj.CON.wheel2.udot.max];
            u_dot_min = [obj.CON.wheel1.udot.min; obj.CON.wheel2.udot.min];
            U_dot_max = repmat(u_dot_max,N-1,1);
            U_dot_min = repmat(u_dot_min,N-1,1);
            u_star = u_trj(:);

            L1 = [D; -D];
            q1 = [U_dot_max - D*u_star;
                 -U_dot_min + D*u_star];
             
            % state constraints
            % form the constraint matrix C
            C = cell(1,N);
            m = [1 0 0; 0 1 0];
            [C{:}] = deal(m);
            C = blkdiag(C{:});
            C = [C; -C];
            
            x_min = obj.CON.state.x.min - s(1,2:end);
            x_max = obj.CON.state.x.max - s(1,2:end);
            y_min = obj.CON.state.y.min - s(2,2:end);
            y_max = obj.CON.state.y.max - s(2,2:end);

            x_con_max = [x_max(:), y_max(:)];
            x_con_max = x_con_max';
            x_con_max = x_con_max(:);
            x_con_min = [x_min(:), y_min(:)];
            x_con_min = x_con_min';
            x_con_min = x_con_min(:);
            x_con = [x_con_max; -x_con_min];
            
            F = ilc.F;
            d = ilc.filter.x;
            L2 = C*F;
            q2 = x_con - C*d; 
            
            % combine input and state constraints
            L = [L1; L2];
            q = [q1; q2];
            
            
        end
        
    end
end