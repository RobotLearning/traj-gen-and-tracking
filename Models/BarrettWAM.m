% Anthropomorphic arm (3shoulder+1elbow+1wrist)

classdef BarrettWAM < Robot

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
        % reference shown in joint space?
        flag_ref_jsp
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            %assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
            
            % set observation matrix
            obj.C = STR.C;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
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
            %assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
            %       'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
            obj.flag_jspace = ~sim.cartesian;
            obj.flag_ref_jsp = sim.jref;
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
        % TODO: divide into several methods?
        function obj = BarrettWAM(par,con,cost,sim)
            
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
        
        % provides nominal model
        function [x_dot,varargout] = nominal(obj,~,x,u,flg)
            % differential equation of the forward dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            [x_dot,dfdx,dfdu] = obj.dynamics(x,u,true);
            
            % return jacobian matrices
            if flg
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            end
        end
        
        % provides actual model
        function x_dot = actual(obj,~,x,u)
            
            % differential equation of the inverse dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            
            % change the masses slightly
            par = obj.PAR;
            %%{
            par.links(1).m = 0.00000 + .5; 
            par.links(2).m = 0.00000 + .5;
            par.links(3).m = 3.53923 + .5; 
            par.links(4).m = 1.03409 + .5;
            par.links(5).m = 2.28843 + .1;
            par.links(6).m = 0.25655 + .1; 
            par.links(7).m = 0.63285 + .1;
            %}
            x_dot = barrettWamDynamicsArt(x,u,par,false);
            
            
        end
        
        % run kinematics using an external function
        function [x1,x2] = kinematics(obj,q)
            
            %TODO:
        end
        
        % call inverse kinematics from outside
        function q = invKinematics(obj,x)
            
            %TODO:
        end
                   
        % dynamics to get u
        function u = invDynamics(obj,q,qd,qdd)
            % inverse dynamics model taken from SL
            u = barrettWamInvDynamicsNE(q,qd,qdd,obj.PAR);
            %u = barrettWamInvDynamicsArt(q,qd,qdd,obj.PAR);
        end
        
        % dynamics to qet Qd = [qd,qdd]
        function [Qd, varargout] = dynamics(obj,Q,u,flag)
            if flag
                [Qd, dfdx, dfdu] = barrettWamDynamicsArt(Q,u,obj.PAR,flag);
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            else
                Qd = barrettWamDynamicsArt(Q,u,obj.PAR,flag);
            end
        end
        
        % make an animation of the robot manipulator
        function animateArm(obj,q_actual,s)
            %TODO:
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            %TODO:
            
        end
        
    end
end