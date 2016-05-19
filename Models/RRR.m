% Two-link planar revolute robot manipulator (RR)

classdef RRR < Robot

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
                         
            % check that the input has all the fields
            % TODO: is there a better way?
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
            %assert(all(strcmp(fieldnames(obj.SIM), fieldnames(sim))));
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            obj.SIM.eps_m = sim.eps_m;
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
        function obj = RRR(par,con,cost,sim)
            
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
        function x_dot = actual(obj,~,x,u)
            
            % assume we dont know the actuator inertia
            par.const.g = obj.PAR.const.g;
            par.link1.motor.inertia = 0.10;%obj.PAR.link1.motor.inertia;
            par.link2.motor.inertia = 0.08;
            par.link3.motor.inertia = 0.15;
            
            % differential equation of the inverse dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            x_dot = RRRDynamics(x,u,par,false);
            
        end
        
        %% Kinematics and dynamics models
        
        % run kinematics using an external function
        function [x1,x2,x3,Ahmats] = kinematics(obj,q,varargin)
            
            [x1,x2,x3,Ahmats] = RRRKinematics(q,obj.PAR);
            
            if nargin == 3
                angle = varargin{1};
                [x1,x2,x3,Ahmats] = obj.rotate(angle,x1,x2,x3,Ahmats);
            end
        end
        
        function [x1,x2,x3,mats] = rotate(obj,phi,x1,x2,x3,mats)
            % rotate everything by phi degrees
            R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
            x1 = R * x1;
            x2 = R * x2;
            x3 = R * x3;
            for timeIdx = 1:size(mats,4)
                Ahmats = mats(:,:,:,timeIdx);
                for slice = 1:size(Ahmats,3)
                    Ahmats(1:2,1:2,slice) = R * Ahmats(1:2,1:2,slice);
                end
                mats(:,:,:,timeIdx) = Ahmats;
            end
        end
                   
        % dynamics to get u
        % TODO
        function u = invDynamics(obj,q,qd,qdd)
            u = RRRInvDynamics(q,qd,qdd,obj.PAR);
        end
        
        % dynamics to qet Qd = [qd,qdd]
        % TODO
        function [Qd, varargout] = dynamics(obj,Q,u,flag)
            if flag
                [Qd, dfdx, dfdu] = RRRDynamics(Q,u,obj.PAR,flag);
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            else
                Qd = RRRDynamics(Q,u,obj.PAR,flag);
            end
        end
        
        % get end effector position
        function [x,xd,mats] = getEndEffectorState(obj,q,qd,varargin)
            
            dofs = size(q,1);
            lenq = size(q,2);
            if nargin == 4
                phi = varargin{:};
                R = [cos(phi) -sin(phi); 
                     sin(phi) cos(phi)];
            else
                R = eye(2);
            end      
            
            [~,~,x,mats] = obj.kinematics(q,varargin{:});            
            xd = zeros(2,lenq);
            for i = 1:lenq
                obj.calcJacobian(q(:,i));
                xd(:,i) = R*obj.jac*qd(:,i);
            end
        end
        
        % get jacobian at current q
        function calcJacobian(obj,q)
            
            l1 = obj.PAR.link1.length;
            l2 = obj.PAR.link2.length;
            l3 = obj.PAR.link3.length;
            obj.jac = [-l1*sin(q(1)) - l2*sin(q(1)+q(2)) - l3*sin(q(1)+q(2)+q(3)), ...
                       -l2*sin(q(1)+q(2)) - l3*sin(q(1)+q(2)+q(3)), ...
                       -l3*sin(q(1)+q(2)+q(3)); 
                       l1*cos(q(1)) + l2*cos(q(1)+q(2)) + l3*cos(q(1)+q(2)+q(3)), ...
                       l2*cos(q(1)+q(2)) + l3*cos(q(1)+q(2)+q(3)), ...
                       l3*cos(q(1)+q(2)+q(3))];
        end
        
        %% Inverse Kinematics functions
        
        % call inverse kinematics from outside
        function q = invKinematics(obj,x,var)
            
            if isscalar(var)
                % phi must have been provided
                phi = var;
                q = RRRInvKinematics(x,phi,obj.PAR);
            else
                Ahmat = var;
                % Ahmat is the last hom. matrix for the last joint
                % get the approach vector
                a = Ahmat(1:3,1);
                % get phi
                phi = Atan2(a(2),a(1));
                q = RRRInvKinematics(x,phi,obj.PAR);
            end
        end        
        
        % inverse kin for table tennis
        function [qf,qfdot] = invKinTableTennis(obj,Q0,racket)
            
            dof = length(Q0)/2;
            q0 = Q0(1:dof);
            % rotate some variables for drawing in 2D simulation
            R = [0 -1; 1 0];
            % feed to inverse kinematics to get qf
            try
                normalRot = R*racket.normal;
                phiVHP = atan2(normalRot(2),normalRot(1));
                qf = obj.invKinematics(R*racket.pos, phiVHP);
                obj.calcJacobian(qf);
                qfdot = obj.jac \ (R*racket.vel);
            catch ME
                disp('Virtual Hitting Point outside of workspace');
                qf = q0;
                qfdot = zeros(dof,1);
            end
        end
        
        %% Check safety here
        
        function [q,qd,qdd] = checkJointLimits(obj,q,qd,qdd)
            
            dofs = size(q,1);
            len = size(q,2);
            con = obj.CON;
            umax = repmat(con.u.max,1,len);
            umin = repmat(con.u.min,1,len);
            qmax = repmat(con.q.max,1,len);            
            qmin = repmat(con.q.min,1,len);
            qdmax = repmat(con.qd.max,1,len);
            qdmin = repmat(con.qd.min,1,len);
            qddmax = repmat(con.qdd.max,1,len);
            qddmin = repmat(con.qdd.min,1,len); 
            u = zeros(dofs,len);
%             for i = 1:len
%                 u(:,i) = obj.invDynamics(q(:,i),qd(:,i),qdd(:,i));
%             end
            
            %obj.displayMaxInfo(q,qd,qdd,u);
            
            try              
                assert(~any(any((q > qmax | q < qmin))),'Joint limits violated!');
                assert(~any(any((qd > qdmax | qd < qdmin))), 'Vel limits violated!');
                assert(~any(any((qdd > qddmax | qdd < qddmin))), 'Acc limits violated!');
                assert(~any(any((u > umax | u < umin))), 'Torque limits violated!');
            catch ME
                disp(ME.message);
                disp('Not moving the robot!');
                q0 = q(:,end);
                q = q0 * ones(1,len);
                qd = zeros(dofs,len);
                qdd = zeros(dofs,len);
            end
        end
        
        % Display max and min info about trajectory and control inputs
        function displayMaxInfo(obj,q,qd,qdd,u)
            
            dofs = size(q,1);
            qmax = max(q,[],2);
            qmin = min(q,[],2);
            qdmax = max(qd,[],2);
            qdmin = min(qd,[],2);
            qddmax = max(qdd,[],2);
            qddmin = min(qdd,[],2);
            umax = max(u,[],2);
            umin = min(u,[],2);
            
            fprintf('qmax 	 qdmax 	 qddmax    umax\n');
            for i = 1:dofs
                fprintf('%-8.2f %-8.2f %-8.2f %-8.2f\n',...
                    qmax(i),qdmax(i),qddmax(i),umax(i));
            end
            fprintf('qmin 	 qdmin 	 qddmin    umin\n');
            for i = 1:dofs
                fprintf('%-8.2f %-8.2f %-8.2f %-8.2f\n',...
                    qmin(i),qdmin(i),qddmin(i),umin(i));
            end
                
        end
        
        function q = clampJointLimits(obj,q)
            
            con = obj.CON;
            q = min(con.q.max,q);
            q = max(con.q.min,q);
        end
        
        function q = checkContactWithTable(obj,q)
            % TODO:
        end        
        
        %% Drawing functions here
        % to draw the robots joints and endeffector 
        % for one posture only
        function [joints,endeff,racket] = drawPosture(obj,q,varargin)
            
            if nargin == 3
                rotAngle = varargin{:};
            end
            
            [x1,x2,x3,mats] = obj.kinematics(q,rotAngle);
            
            base = [0,0];
            joints = [base; x1'; x2'];
            endeff = x3';
            
            racket_orient = mats(1:2,1,3);
            racket_dir = mats(1:2,2,3);
            racket_radius = 0.08;
            racket1 = x3 - racket_radius * racket_dir;
            racket2 = x3 + racket_radius * racket_dir;
            racket = [racket1'; racket2'];
            
        end
        
        % Make an animation of the robot manipulator
        function animateArm(obj,q_actual,s)
            [x1,x2,x3,~] = obj.kinematics(q_actual);
            if obj.flag_ref_jsp
                [~,s] = obj.kinematics(s);
            end
            animateRRR(x1,x2,x3,s);
        end
        
        
        %% get lifted model constraints
        % now only aILC may use it
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            h = obj.SIM.h;
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            umin(1,:) = obj.CON.link1.u.min - u_trj(1,:);
            umin(2,:) = obj.CON.link2.u.min - u_trj(2,:);
            umin(3,:) = obj.CON.link3.u.min - u_trj(2,:);
            umax(1,:) = obj.CON.link1.u.max - u_trj(1,:);
            umax(2,:) = obj.CON.link2.u.max - u_trj(2,:);
            umax(3,:) = obj.CON.link3.u.max - u_trj(3,:);

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/h;
            D = D(1:end-dimu,:);
            
            u_dot_max = [obj.CON.link1.udot.max; 
                         obj.CON.link2.udot.max;
                         obj.CON.link3.udot.max];
            u_dot_min = [obj.CON.link1.udot.min; 
                         obj.CON.link2.udot.min;
                         obj.CON.link3.udot.min];
            U_dot_max = repmat(u_dot_max,N-1,1);
            U_dot_min = repmat(u_dot_min,N-1,1);
            u_star = u_trj(:);

            L = [D; -D];
            q = [U_dot_max - D*u_star;
                 -U_dot_min + D*u_star];
            
        end
        
    end
end