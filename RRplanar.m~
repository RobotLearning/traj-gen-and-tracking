% Two-link planar revolute robot manipulator (RR)

classdef RRplanar < Robot

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
        % joint space coordinates of the nominal trajectory
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
            obj.CON.link1.qd.max = Inf;
            obj.CON.link1.qd.min = -Inf;
            obj.CON.link1.qdd.max = Inf;
            obj.CON.link1.qdd.min = -Inf;
            obj.CON.link1.u.max = Inf;
            obj.CON.link1.u.min = -Inf;
            obj.CON.link2.q.max = Inf;
            obj.CON.link2.q.min = -Inf;
            obj.CON.link2.qd.max = Inf;
            obj.CON.link2.qd.min = -Inf;
            obj.CON.link2.qdd.max = Inf;
            obj.CON.link2.qdd.min = -Inf;            
            obj.CON.link2.u.max = Inf;
            obj.CON.link2.u.min = -Inf;
            
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.dimx = 2;
            obj.SIM.dimu = 2;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
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
        function obj = RRplanar(par,con,cost,sim)
            
            obj.SIM = sim;            
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;        
            % cost function handle
            obj.COST = cost.Q;
        end
        
        % provides nominal model
        function [x_dot,varargout] = nominal(~,obj,x,u,flg)
            % differential equation of the inverse dynamics
            % x_dot = Ax + B(x)u + C
            [x_dot,dfdx,dfdu] = obj.inverseDynamics(x,u,true);
            
            % return jacobian matrices
            if flg
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            end
        end
        
        % provides actual model
        function x_dot = actual(~,obj,x,u)
            
            % change parameters
            par.const.g = obj.PAR.const.g;
            par.link1.mass = 1.0 * obj.PAR.link1.mass;
            par.link2.mass = 1.0 * obj.PAR.link2.mass;
            par.link1.length = 1.2 * obj.PAR.link1.length;
            par.link2.length = 1.2 * obj.PAR.link2.length;
            par.link1.centre.dist = obj.PAR.link1.centre.dist;
            par.link2.centre.dist = obj.PAR.link2.centre.dist;
            par.link1.inertia = 1.0 * obj.PAR.link1.inertia;
            par.link2.inertia = 1.0 * obj.PAR.link2.inertia;
            par.link1.motor.inertia = 1.0 * obj.PAR.link1.motor.inertia;
            par.link2.motor.inertia = 1.0 * obj.PAR.link2.motor.inertia;
            par.link1.motor.gear_ratio = 1.0 * obj.PAR.link1.motor.gear_ratio;
            par.link2.motor.gear_ratio = 1.0 * obj.PAR.link2.motor.gear_ratio;
            obj.PAR = par;
            
            % differential equation of the inverse dynamics
            % x_dot = Ax + B(x)u + C
            x_dot = obj.inverseDynamics(x,u,false);
            
        end
        
        % run kinematics using an external function
        function [x1,x2] = kinematics(obj,q)
            
            [x1,x2] = RRplanarKinematics(q,obj.PAR);
        end
        
        % call inverse kinematics from outside
        function q = inverseKinematics(obj,x)
            
            q = RRplanarInverseKinematics(x,obj.PAR);
        end
                   
        % dynamics to get u
        function u = dynamics(obj,q,qd,qdd)
            u = RRplanarDynamics(q,qd,qdd,obj.PAR);
        end
        
        % inverse dynamics to qet Qd = [qd,qdd]
        function Qd = inverseDynamics(obj,Q,u,flag)
            Qd = RRplanarInverseDynamics(Q,u,obj.PAR,flag);
        end
        
        % make an animation of the robot manipulator
        function animateArm(obj,q_actual,s)
            [x1,x2] = obj.kinematics(q_actual);
            animateRR(x1,x2,s);
        end

        % using a simple inverse kinematics method
        % TODO: extend using planning to incorporate constraints
        function Traj = trajectory(obj,t,x_des)

            h = obj.SIM.h;
            obj.q = RRplanarInverseKinematics(x_des,obj.PAR);
            % check for correctness
            %[x1,x2] = RRplanarKinematics(q,PAR);
            obj.qd = diff(obj.q')' / h; 
            obj.qdd = diff(obj.qd')' / h; 

            % keep velocity and acceleration vectors the same length as displacements
            obj.qd(:,end+1) = obj.qd(:,end);
            obj.qdd(:,end+1) = obj.qdd(:,end);
            obj.qdd(:,end+1) = obj.qdd(:,end);

            % get the desired inputs
            ud = zeros(size(obj.q));
            for i = 1:size(obj.q,2)
                ud(:,i) = RRplanarDynamics(obj.q(:,i),obj.qd(:,i),...
                                           obj.qdd(:,i),obj.PAR);
            end
            
            Traj = Trajectory(t,[],x_des,ud);
        end
        
    end
end