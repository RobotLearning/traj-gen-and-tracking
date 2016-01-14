% Robot abstract class is the parent for particular
% robot arm manipulators.

classdef (Abstract) Robot < Model
        
    % Field necessary for all robotic manipulators?
    properties (Abstract)        

        % jacobian
        jac
        % flag for representing learned trajectories
        % in operational (0) or joint space (1) ?
        flag_jspace
        % flag indicating joint space references
        flag_ref_jsp
    end
    
    % methods to be implemented
    methods (Abstract)
        
        % kinematics
        kinematics(q)
        % inverse kinematics
        invKinematics(x)
        % inverse dynamics to get u
        invDynamics(q,qd,qdd)
        % (direct) dynamics to qet Qd = [qd,qdd]
        dynamics(Q,u)
        % make an animation of the robot manipulator
        animateArm(x)
        
    end
    
    % methods that robots share
    methods (Access = public)
        
        %% Generate table tennis trajectories with the VHP method
        % TODO
        function Traj = generateTTTwithVHP(obj,VHP,ballPred,ballTime,ballDes,time2reach)
            
            ballPredWithTime = [ballPred;ballTime];
            ballAtVHP = interp1(ballPredWithTime(1,:)',ballPredWithTime(2:5,:)',VHP);
            timeAtVHP = ballAtVHP(end);
            ballAtVHP = [VHP;ballAtVHP(1:end-1)'];
            ballPosAtVHP = ballAtVHP(1:2);
            ballInVelAtVHP = ballAtVHP(3:4);
            
            % GET DESIRED OUTGOING VELOCITY OF THE BALL AT VHP
            
            loadTennisTableValues();
            % approximate with a linear model
            ballOutVelAtVHP(1) = (ballDes(1) - ballPosAtVHP(1))/time2reach;
            ballOutVelAtVHP(2) = (ballDes(2) - ballPosAtVHP(2) - ...
                                0.5*gravity*time2reach^2)/time2reach;
            
            % initialize using a linear model (no drag)
            linFlightTraj = @(t) [ballPosAtVHP + ballOutVelAtVHP(:)*t;
                                  ballOutVelAtVHP(:)] + ...
                                 [0;0.5*gravity*t^2;0;gravity*t];
            flightModel = @(t,x) [x(3);
                                  x(4);                                  
                              -Cdrag * x(3) * sqrt(x(3)^2 + x(4)^2);
                      gravity - Cdrag * x(4) * sqrt(x(3)^2 + x(4)^2)];
            % boundary value condition
            bc = @(x0,xf) [x0(1) - ballPosAtVHP(1);
                           x0(2) - ballPosAtVHP(2);
                           xf(1) - ballDes(1);
                           xf(2) - ballDes(2)];
            meshpoints = 50;
            solinit = bvpinit(linspace(0,time2reach,meshpoints),...
                               linFlightTraj);
            sol = bvp4c(flightModel,bc,solinit);
            ballOut = deval(sol,0);
            ballOutVelAtVHP = ballOut(3:4);
            
            % GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
            
            % inverting the mirror law
            normal = (ballOutVelAtVHP - ballInVelAtVHP) ...
                          ./ norm(ballOutVelAtVHP - ballInVelAtVHP);
            velOutAlongNormal = ballOutVelAtVHP' * normal;
            velInAlongNormal = ballInVelAtVHP' * normal;
            racketVelAlongNormal = (velOutAlongNormal + ... 
                CRR * velInAlongNormal) / (1 + CRR);
            % racket velocity along racket plane we take to be zero
            % this implies no spinning
            racketVel = racketVelAlongNormal * normal;
                            
            %phiVHP = atan2(normal(2),normal(1))-pi/2;
            % get desired orientation angle phi at VHP
            normalRot = R'*normal;
            phiVHP = atan2(normalRot(2),normalRot(1));
            % feed to inverse kinematics to get qf
            try
                qf = RRRInvKinematics(R'*ballPosAtVHP,phiVHP,rrr.PAR);
                rrr.calcJacobian(qf);
                qfdot = rrr.jac \ (R'*racketVel);
            catch ME
                disp('Virtual Hitting Point outside of workspace');
                qf = q0;
                qfdot = zeros(3,1);
            end

            q0dot = zeros(3,1);
            
            % GET 3RD DEGREE POLYNOMIALS
            thit = timeAtVHP;
            t = dt:dt:thit;
            t2 = dt:dt:time2return;
            M = @(t0,tf) [t0^3 t0^2 t0^1 1;
                          3*t0^2 2*t0 1 0;
                          tf^3 tf^2 tf^1 1;
                          3*tf^2 2*tf 1 0];
            qStrike = zeros(3,length(t));
            qReturn = zeros(3,length(t2));            
            for m = 1:3
                %q0dot is zero
                Qstrike = [q0(m); q0dot(m); qf(m); qfdot(m)]; % strike
                Qreturn = [qf(m); qfdot(m); q0(m); q0dot(m)]; % return
                a = M(0,thit) \ Qstrike;
                b = M(0,time2return) \ Qreturn;
                qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
                qReturn(m,:) = b(1)*t2.^3 + b(2)*t2.^2 + b(3)*t2 + b(4);
            end
            q = [qStrike,qReturn];
        end
        
        %% Generate inputs for a reference
        % using a simple inverse kinematics method
        % if flag is 1 then reference in joint space!
        function Traj = generateInputs(obj,t,ref)

            h = obj.SIM.h;
            dt = t(2) - t(1); % different from h if downsampled
            dim = obj.SIM.dimx / 2;
            dimu = obj.SIM.dimu;
            Cout = obj.C;
            
            if (obj.flag_ref_jsp)
                q = ref(1:dim,:);
                qd = ref(dim+1:2*dim,:);
            else
                % assuming that only positions are demonstrated
                q = obj.invKinematics(ref);
                qd = diff(q')' / dt; 
                % start with zero initial velocity
                %qd = [zeros(dim,1), qd];
                %qdd = [zeros(dim,1), qdd];
                % assume you end with same velocity as before
                %qd(:,end+1) = qd(:,end);
                % linearly extrapolate (zero jerk at traj end)
                qd(:,end+1) = qd(:,end) + qd(:,end) - qd(:,end-1);
            end
            
            qdd = diff(qd')' / dt; 
            % assume you end with same acceleration as before
            %qdd(:,end+1) = qdd(:,end);
            % this leads to large decelerations!

            % get the desired inputs
            Nu = length(t) - 1;
            uff = zeros(dimu,Nu);
            for i = 1:Nu
                uff(:,i) = obj.invDynamics(q(:,i),qd(:,i),qdd(:,i));
            end
            
            % check for joint space representation
            if obj.flag_jspace == 1
                Traj = Trajectory(t,Cout*[q;qd],uff,[]);
            else
                %xd = obj.jac * qd;
                x =  ref;
                xd = diff(x')'/dt;
                xd(:,end+1) = xd(:,end);
                Traj = Trajectory(t,Cout*[x;xd],uff,[]);
            end            
        end     
        
        % Method useful when modifying DMPs directly
        % varargin is for different initial conditions
        function [Traj,dmps] = generateInputsWithDMP(obj,t,numbf,ref,varargin)

            h = obj.SIM.h;
            dim = obj.SIM.dimx / 2;
            dimu = obj.SIM.dimu;
            Cout = obj.C;
            
            assert(size(Cout,2)==obj.SIM.dimx,...
                'Currently works only for full observation');
            
            if (obj.flag_ref_jsp)
                qdes = ref(1:dim,:);
                qddes = ref(dim+1:2*dim,:);
            else
                % assuming that only positions are demonstrated
                qdes = obj.invKinematics(ref);
                qddes = diff(qdes')' / h; 
                % assume you end with same velocity as before
                qddes(:,end+1) = qddes(:,end);
            end
            
            % make DMPs that smoothens reference, one for each output
            
            % this is quite unnecessary as velocities are always zero in
            % goal positions
            goal = [qdes(:,end);qddes(:,end)];
            if nargin > 4
                yin = varargin{1};
                yin = reshape(yin,dim,2);
            else
                yin = [qdes(:,1),qddes(:,1)];     
            end
            [dmps,s] = obj.dmpTrajectory(t,numbf,goal,yin,[qdes;qddes]);
            
            q = s(1:dim,:);
            qd = s(dim+1:end,:);            
            qdd = diff(qd')' / h; 
            % assume you end with same acceleration as before
            qdd(:,end+1) = qdd(:,end);

            % get the desired inputs
            Nu = length(t) - 1;
            uff = zeros(dimu,Nu);
            for i = 1:Nu
                uff(:,i) = obj.invDynamics(q(:,i),qd(:,i),qdd(:,i));
            end
            
            % check for joint space representation
            if obj.flag_jspace == 1
                Traj = Trajectory(t,Cout*[q;qd],uff,[]);
            else
                %xd = obj.jac * qd;
                x =  ref;
                xd = diff(x')'/h;
                xd(:,end+1) = xd(:,end);
                Traj = Trajectory(t,Cout*[x;xd],uff,[]);
            end            
        end
        
        % Create inputs and trajectory class from a DMP
        % DMP is assumed to be in joint space
        function traj = generateInputsForDMP(obj,dmp,N)

            h = obj.SIM.h;
            dim = obj.SIM.dimx / 2;
            dimu = obj.SIM.dimu;
            Cout = obj.C;
            t = h * (1:N);
            assert(dmp(1).can.dt == obj.SIM.h, 'Time steps should be equal!');
            assert(size(Cout,2)==obj.SIM.dimx,...
                'Currently works only for full observation');
            
            q = zeros(dim,N);
            qd = zeros(dim,N);
            for i = 1:length(dmp)
                % make sure they are reset
                dmp(i).resetStates();
                [~,si] = dmp(i).evolve(N);
                q(i,:) = si(1,:);
                qd(i,:) = si(2,:);
            end
       
            qdd = diff(qd')' / h; 
            % assume you end with same acceleration as before
            qdd(:,end+1) = qdd(:,end);

            % get the desired inputs
            Nu = length(t) - 1;
            uff = zeros(dimu,Nu);
            for i = 1:Nu
                uff(:,i) = obj.invDynamics(q(:,i),qd(:,i),qdd(:,i));
            end
            
            % check for joint space representation
            if obj.flag_jspace == 1
                traj = Trajectory(t,Cout*[q;qd],uff,[]);
            else
                %xd = obj.jac * qd;
                x =  ref;
                xd = diff(x')'/h;
                xd(:,end+1) = xd(:,end);
                traj = Trajectory(t,Cout*[x;xd],uff,[]);
            end            
        end
        
        %% generating feedback with LQR to stabilize the robot
        function generateFeedback(obj,traj)
            
            % calculate the optimal feedback law
            t = traj.t;
            Q = obj.COST.Q;
            R = obj.COST.R;
            Cout = obj.C;
            N = length(t)-1;
            h = t(2) - t(1);
            % get linear time variant matrices around trajectory
            [Ad,Bd] = obj.linearize(traj);
            lqr = LQR(Q,R,Q,Ad,Bd,Cout,N,h,true);
            K = lqr.computeFinHorizonLTV();
            
            traj.K = K;
            
        end
        
    end
        
    
end