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
        function obj = BarrettWAM(par,con,cost,sim)
            
            obj.SIM = sim;            
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;        
            % cost function handle
            obj.COST = cost;
            % jacobian in case we need it
            obj.jac = [];
        end
        
        %% Dynamics functions here
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
        function x_dot = actual(obj,t,x,u)
            
            % differential equation of the inverse dynamics
            % x_dot = A(x) + B(x)u
            
            % load actual values
            loadActualBarrettValues;
            par = obj.PAR;
            par.link0 = link0;
            par.links = links;
            d = @(t) 0.0; %0.2 + 0.5 * sin(10*t);
            x_dot = barrettWamDynamicsArt(x,u,par,false) + d(t);            
            
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
        
        %% Kinematics related functions here
        % Calculate the racket orientation based on quaternion
        function racketOrient = calcRacketOrientation(obj,cartOrient)
            
            % quaternion transformation of pi/4 from endeff to racket
            rot = [cos(pi/4); -sin(pi/4); 0; 0];
            
            racketOrient = mult2Quat(cartOrient,rot);
            
%             racketOrient(1) = cartOrient(1)*rot(1) - ...
%                               cartOrient(2:4)*rot(2:4);
%             racketOrient(2) = cartOrient(2)*rot(1) + ...
%                               cartOrient(1)*rot(2) + ...
%                               cartOrient(3)*rot(4) - ...
%                               cartOrient(4)*rot(3);
%             racketOrient(3) = cartOrient(1)*rot(3) - ...
%                               cartOrient(2)*rot(4) + ...
%                               cartOrient(3)*rot(1) + ...
%                               cartOrient(4)*rot(2);
%             racketOrient(4) = cartOrient(1)*rot(4) + ...
%                               cartOrient(2)*rot(3) - ...
%                               cartOrient(3)*rot(2) + ...
%                               cartOrient(4)*rot(1);
            
        end
        
        % run kinematics using an external function
        % return endeffector coordinates, vel and orientation
        % TODO: should return what barrettWAMKinematics returns
        function [x,xd,o] = kinematics(obj,Q)
              
            assert(size(Q,1) == 14, 'velocities not fed in!');
            dim = size(Q,1)/2;
            lenq = size(Q,2);
            q = Q(1:dim,:);
            qd = Q(dim+1:end,:);
            x = zeros(3,lenq);
            xd = zeros(3,lenq);
            o = zeros(4,lenq);
            for i = 1:lenq
                [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q(:,i),obj.PAR);
                o(:,i) = rot2Quat(Amats(6,1:3,1:3));
                obj.jac = jacobian(xLink,xOrigin,xAxis);
                x(:,i) = xLink(6,:)';
                xd(:,i) = obj.jac(1:3,:) * qd(:,i);
            end
        end
        
        % racket configurations and velocities are returned
        % racket orientations are also returned as a quaternions        
        function [x,xd,o] = calcRacketState(obj,Q)
            
            assert(size(Q,1) == 14, 'velocities not fed in!');
            dim = size(Q,1)/2;
            lenq = size(Q,2);
            q = Q(1:dim,:);
            qd = Q(dim+1:end,:);
            x = zeros(3,lenq);
            xd = zeros(3,lenq);
            o = zeros(4,lenq);
            for i = 1:lenq
                [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q(:,i),obj.PAR);
                quat = rot2Quat(Amats(6,1:3,1:3));
                o(:,i) = obj.calcRacketOrientation(quat);
                obj.jac = jacobian(xLink,xOrigin,xAxis);
                x(:,i) = xLink(6,:)';
                xd(:,i) = obj.jac(1:3,:) * qd(:,i);
            end
        end
        
        % calculate the geometric jacobian using the common jacabian function
        function jac = calcJacobian(obj,q)
            
            [xLink,xOrigin,xAxis,~] = barrettWamKinematics(q,obj.PAR);
            jac = jacobian(xLink,xOrigin,xAxis);
            obj.jac = jac; 
        end
        
        %% Inverse Kinematics functions
        
        % Inverse Kinematics for table tennis
        % Specifically, trying to keep initial slide of the racket
        function [qf,qfdot] = invKinTableTennis1(obj,Q0,racket)
            
            racketPos = racket.pos;
            racketNormal = racket.normal;
            racketVel = racket.vel;
            racketAngularVel = racket.angvel;
            
            try
                tic;
                % get the slide of the original racket orientation
                [~,~,o] = obj.calcRacketState(Q0);
                rotMatrix0 = quat2Rot(o);
                slide0 = rotMatrix0(1:3,2);    
                slide0 = slide0./norm(slide0,2);
                % add a comfortable slide close to original slide
                a = racketNormal;
                % project slide0 to the racket plane
                %projNormal = racketNormal*racketNormal';
                %projRacket = eye(3) - projNormal;
                %s = projRacket * slide0;    
                s = slide0 - (racketNormal'*slide0).*racketNormal;    
                assert(abs(norm(s,2)-1) < 1e-2,'slide not normalized');
                s = s./norm(s,2);
                assert(s'*a < 1e-3,'slide calculation not working!');
                n = crossProd(s,a);                
                rotMatrix = [n,s,a];
                quatRacket = rot2Quat(rotMatrix);
                % get endeffector position and quaternion
                ePos = racketPos;
                rotBack = [cos(-pi/4); -sin(-pi/4); 0; 0];
                eQuat = mult2Quat(quatRacket,rotBack);
                q0 = Q0(1:7);
                qf = obj.invKinematics(ePos(:),eQuat(:),q0(:));
                timeInvKin = toc;
                fprintf('InvKin took %f sec.\n',timeInvKin);
                obj.calcJacobian(qf);                
                qfdot = obj.jac \ [racketVel;racketAngularVel];
            catch ME
                disp(ME.message);
                disp('InvKin problem. Not moving the robot...');                
                %disp('Virtual Hitting Point outside of workspace');
                qf = q0;
                qfdot = zeros(dof,1);
            end
        end
        
        % call inverse kinematics from outside
        % o is a quaternion
        function q = invKinematics(obj,x,o,q0)
            
            % not masking any cartesian space variables
            %M = [1 1 1 1 1 1];
    
            %fprintf('Computing IK for WAM...\n'); tic;
            for k = 1:1:size(x,2) 
        
                T = [quat2Rot(o(:,k)),x(:,k);zeros(1,3),1];
                q(:,k) = invKinematics(T, q0, obj.PAR);
                q0 = q(:,k); % use previous joint values as initial guess for ikine
            end
            %fprintf('finished IK in %g seconds.\n', toc);
            
            % wrap to [-pi,pi];
            %idx = (q < -pi) | (pi < q);
            %lhs = q > -pi; % left hand side limit reached
            %q(idx) = mod(q(idx) + pi, 2*pi) - pi; %mod to [-pi,pi]
            %q(idx & (q == -pi) & lhs) = pi;
            
            %q = unwrap(q);
            
        end        

        %% Check safety here
        
        function checkJointLimits(obj,q,qd,qdd)
            
            len = size(q,2);
            con = obj.CON;
            qmax = repmat(con.q.max,1,len);            
            qmin = repmat(con.q.min,1,len);
            qdmax = repmat(con.qd.max,1,len);
            qdmin = repmat(con.qd.min,1,len);
            qddmax = repmat(con.qdd.max,1,len);
            qddmin = repmat(con.qdd.min,1,len);            
            assert(sum(sum((q > qmax) + (q < qmin))) == 0,'joint limits violated');
            assert(sum(sum((qd > qdmax) + (qd < qdmin))) == 0, 'vel limits violated');
            assert(sum(sum((qdd > qddmax) + (qdd < qddmin))) == 0, 'acc limits violated');
             
        end
        
        %% Drawing functions here
        % to draw the robots joints and endeffector 
        % for one posture only
        function [joints,endeff,racket] = drawPosture(obj,q)
            
            [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q,obj.PAR);
            quat = rot2Quat(Amats(6,1:3,1:3));
            orient = obj.calcRacketOrientation(quat);
            R = quat2Rot(orient(:));
            joints = xOrigin;
            endeff = xLink(6,:);
            
            %x and y are the coordinates of the center of the circle
            %r is the radius of the circle
            %0.01 is the angle step, bigger values will draw the circle faster but
            %you might notice imperfections (not very smooth)
            racket_radius = 0.076;
            ang = 0:0.01:2*pi; 
            racket_x = racket_radius * cos(ang);
            racket_y = racket_radius * sin(ang);
            % transform into base coord.
            racket = repmat(endeff(:),1,length(ang)) + ...
                R * [racket_x; racket_y; zeros(1,length(ang))];
            
        end
        
        % make an animation of the robot manipulator
        %TODO
        function animateArm(obj,qs)
            
            dim = 7;
            assert(size(qs,1) == dim, 'velocities not necessary!');
            NCART = 3;
            x = zeros(NCART,lenq);
            o = zeros(4,lenq);
            for i = 1:lenq
                [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(qs(:,i),obj.PAR);
                quat = rot2Quat(Amats(6,1:3,1:3));
                o(:,i) = obj.calcRacketOrientation(quat);
                x(:,:,i) = [xOrigin;xLink(6,:)];
            end
            
            animateWAM(x,o);
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            %TODO:
            
        end
        
    end
end