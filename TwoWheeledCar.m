% Two-wheeled 2D car or robot kinematical model 
% TODO: extend for dynamics

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
        % cartesian coordinates of the nominal trajectory
        x, xd, xdd
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            % initialize everything to zero
            obj.PAR.wheel1.radius = 0;
            obj.PAR.wheel2.radius = 0;
            obj.PAR.length = 0;
                         
            % check that the input has all the fields
            % TODO: is there a better way?
            assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % initialize all fields to zero
            obj.CON.state.x.max = 0;
            obj.CON.state.x.min = 0;
            obj.CON.state.y.max = 0;
            obj.CON.state.y.min = 0;
            obj.CON.state.theta.max = 0;
            obj.CON.state.theta.min = 0;
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
            obj.SIM.dimx = 3;
            obj.SIM.dimu = 2;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        function set.COST(obj, Q)
            obj.COST.Q = Q;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*Q*(x1-x2));
            assert(length(Q) == obj.SIM.dimx);
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
            obj.COST = cost.Q;
        end
        
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
            par.length = obj.PAR.length;
            
            % differential equation of the inverse dynamics
            % x_dot = B(x)u
            x_dot = robotTwoWheelsKinematics(t,x,u,par,false);
            
        end
                           
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
        
        % make an animation of the robot manipulator
        function animate(obj,x,s)
            
            animateCar(x,s,obj.PAR);
        end

        % using a simple inverse kinematics method
        % TODO: extend using planning to incorporate constraints
        function Traj = trajectory(obj,t,x_des)

            N = length(t);
            dimu = obj.SIM.dimu;
            h = obj.SIM.h;
            unom = zeros(dimu,N);
            R1 = obj.PAR.wheel1.radius;
            R2 = obj.PAR.wheel2.radius;
            d = obj.PAR.length;

            for i = 1:N-1
                % discretize B(phi) around phi0
                B = h * [R1/2 * cos(x_des(3,i)), R2/2 * cos(x_des(3,i));
                         R1/2 * sin(x_des(3,i)), R2/2 * sin(x_des(3,i));
                         R1/d,                   -R2/d];
                % deviation along desired trajectory
                delta = x_des(2:3,i+1) - x_des(2:3,i);
                % solve for u
                unom(:,i) = B(2:3,:)\delta;
            end

            unom(:,end) = unom(:,end-1);
            
            Traj = Trajectory(t,[],x_des,unom);
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj)
            
            h = obj.SIM.h;
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
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

            L = [D; -D];
            q = [U_dot_max - D*u_star;
                 -U_dot_min + D*u_star];
            
        end
        
    end
end