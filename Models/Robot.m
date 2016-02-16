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
        
        %% Table Tennis related functions
        
        % Generate optimal tt trajectories
        function [q,qd,qdd] = generateOptimalTTT(obj,racket,ballPred,ballTime,q0)
                  
            dof = length(q0);
            dt = ballTime(2) - ballTime(1);
            time2return = 1.0;
            
            % FOR TESTING SL CONNECTION
            %qest = [2.25;-0.38;-1.27;1.33;-1.86;-0.19;0.77];
            %qdest = [1.10;-0.53;-3.21;1.26;-0.44;-0.09;0.78];
            %timeEst = 0.8;

            %qf = qest;
            %qfdot = qdest;
            %time2hit = timeEst;
            [qf,qfdot,time2hit] = calcOptimalPoly(obj,racket,ballTime,ballPred,q0);
            % round time2hit to nearest dt
            time2hit = dt * ceil(time2hit/dt);
            
            q0dot = zeros(dof,1);
            Q0 = [q0;q0dot];  
            Qf = [qf;qfdot];
            
            % GET 3RD DEGREE POLYNOMIALS            
            pStrike = generatePoly3rd(Q0,Qf,dt,time2hit);
            qStrike = pStrike(1:dof,:);
            qdStrike = pStrike(dof+1:2*dof,:);
            qddStrike = pStrike(2*dof+1:end,:);
            
            pReturn = generatePoly3rd(Qf,Q0,dt,time2return);
            qReturn = pReturn(1:dof,:);
            qdReturn = pReturn(dof+1:2*dof,:);
            qddReturn = pReturn(2*dof+1:end,:);
            
            q = [qStrike,qReturn];
            qd = [qdStrike,qdReturn];
            qdd = [qddStrike,qddReturn];            
            
        end
        
        % Generate 3D table tennis trajectories with the VHP method
        function [q,qd,qdd] = generate3DTTTwithVHP(obj,VHP,ballPred,ballTime,q0)
                        
            loadTennisTableValues();
            dof = length(q0);
            time2reach = 0.5; % time to reach desired point on opponents court
            % land the ball on the centre of opponents court
            ballDes(1) = 0.0;
            ballDes(2) = dist_to_table - 3*table_y/2;
            ballDes(3) = table_z + ball_radius;
            [qf,qfdot,timeAtVHP] = calcPolyAtVHP(obj,VHP,time2reach,ballDes,ballPred,ballTime,q0);

            q0dot = zeros(dof,1);
            Q0 = [q0;q0dot];  
            Qf = [qf;qfdot];            
            dt = ballTime(2) - ballTime(1);
            time2return = 0.5; % time to return to initial configuration
            
            % not moving in case calculation gone wrong
            eps = 1e-2;
            if norm(qf-q0, 2) < eps
                disp('Not moving!');
                totalTime = timeAtVHP + time2return;
                N = floor(totalTime/dt);
                q = q0 * ones(1,N);
                qd = zeros(dof,N);
                qdd = zeros(dof,N);
                return;
            end
            
            % GET 3RD DEGREE POLYNOMIALS            
            pStrike = generatePoly3rd(Q0,Qf,dt,timeAtVHP);
            qStrike = pStrike(1:dof,:);
            qdStrike = pStrike(dof+1:2*dof,:);
            qddStrike = pStrike(2*dof+1:end,:);
            
            pReturn = generatePoly3rd(Qf,Q0,dt,time2return);
            qReturn = pReturn(1:dof,:);
            qdReturn = pReturn(dof+1:2*dof,:);
            qddReturn = pReturn(2*dof+1:end,:);
            
            q = [qStrike,qReturn];
            qd = [qdStrike,qdReturn];
            qdd = [qddStrike,qddReturn];         
            
            % for debugging
            %[x,xd,o] = obj.calcRacketState([q;qd]);
            %rotMs = quat2Rot(o);
            %normals = rotMs(1:3,3,:);
        end
        
        % Generate 2D table tennis trajectories with the VHP method 
        function [q,qd,qdd] = generate2DTTTwithVHP(obj,ballPred,ballTime,q0)
                        
            loadTennisTableValues();
            dof = length(q0);
            
            % define virtual hitting plane (VHP)
            VHP = -0.5;
            time2reach = 0.5; % time to reach desired point on opponents court
            time2return = 0.5; % time to return to initial configuration          
            dt = ballTime(2)-ballTime(1);
            ballFull = [ballPred;ballTime];
            
            % rotate some variables for drawing in 2D simulation
            R = [0 1; -1 0];
            % land the ball on the centre of opponents court
            ballDes(1) = dist_to_table - 3*table_y/2;
            ballDes(2) = table_z + ball_radius;            
            fprintf('Desired landing point: %f\n',ballDes(1));
            
            ballAtVHP = interp1(ballFull(1,:)',ballFull(2:5,:)',VHP);
            timeAtVHP = ballAtVHP(end);
            ballAtVHP = [VHP;ballAtVHP(1:end-1)'];
            ballPosAtVHP = ballAtVHP(1:2);
            ballInVelAtVHP = ballAtVHP(3:4); 
            
            % GET DESIRED OUTGOING VELOCITY OF THE BALL AT VHP            
            ballOutVelAtVHP = calcBallVelOut2D(ballDes,ballPosAtVHP,time2reach);            
            
            % GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
            [racketPos,racketVel,racketOrient] = calcDesRacketState ...
                           (ballPosAtVHP,ballOutVelAtVHP,ballInVelAtVHP);
            
            % feed to inverse kinematics to get qf
            try
                normalRot = R'*racketOrient;
                phiVHP = atan2(normalRot(2),normalRot(1));
                qf = obj.invKinematics(R'*racketPos,phiVHP);
                obj.calcJacobian(qf);
                qfdot = obj.jac \ (R'*racketVel);
            catch ME
                disp('Virtual Hitting Point outside of workspace');
                qf = q0;
                qfdot = zeros(dof,1);
            end
            
            q0dot = zeros(dof,1);
            Q0 = [q0;q0dot];
            Qf = [qf;qfdot];
            
            % GET 3RD DEGREE POLYNOMIALS            
            pStrike = generatePoly3rd(Q0,Qf,dt,timeAtVHP);
            qStrike = pStrike(1:dof,:);
            qdStrike = pStrike(dof+1:2*dof,:);
            qddStrike = pStrike(2*dof+1:end,:);
            
            pReturn = generatePoly3rd(Qf,Q0,dt,time2return);
            qReturn = pReturn(1:dof,:);
            qdReturn = pReturn(dof+1:2*dof,:);
            qddReturn = pReturn(2*dof+1:end,:);
            
            q = [qStrike,qReturn];
            qd = [qdStrike,qdReturn];
            qdd = [qddStrike,qddReturn];
              
        end
        
        %% Generate trajectories and inputs
        
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