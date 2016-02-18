%% Table Tennis class for simulating a match between two robots!

% NOT USED FOR NOW!
classdef TableTennis < handle
    
    properties
        
        % table related parameters
        TABLE
        % net is useful for strategies
        NET
        % ball class
        ball
        % robot 1
        robot1
        % robot 2 (opponent)
        robot2
        % filter for robot 1
        filter1
        % index for robot 1 (useful for animating fast)
        idx
        % initial joint for robot 1
        q0
        % handle structure for drawing animation
        handle
        % finite state automaton for robot strategies
        FSA
        
    end
    
    methods
        
        %% CONSTRUCTOR
        function obj = TableTennis(wam,wam2)
            
            % initialize two robots
            obj.robot1 = wam;
            obj.robot2 = wam2;
            
            % Initialize arm posture
            % initialize the arm with zero velocity on the right hand side
            q0 = [1.8; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
            obj.q0 = q0;
            
            % initialize a ball
            obj.ball = Ball();             
            % initialize animation
            obj.initAnimation(q0);
            % initialize filters
            obj.initFilter();            
            
            obj.FSA = 0; % WAIT
            obj.idx = 1;
            
            loadTennisTableValues();
            % table related values
            obj.TABLE.Z = table_z;
            obj.TABLE.WIDTH = table_width;
            obj.TABLE.LENGTH = table_length;
            % coeff of restitution-friction vector
            obj.TABLE.K = [CFTX; CFTY; -CRT];
            
            % net y value and coeff of restitution of net
            obj.NET.Y = dist_to_table - table_y;
            obj.NET.Xmax = table_width/2 + net_overhang;
            obj.NET.Zmax = table_z + net_height;
            obj.NET.CRN = net_restitution;
            
        end
        
        %% MAIN LOOP
        
        % a robot serves for 2 points
        % and then they switch
        % till one gets 11 points or more with 2 difference
        function match(obj)
            
            tennisBall = obj.ball;
            r1 = obj.robot1;
            r2 = obj.robot2;
            % points for the robots
            points1 = 0;
            points2 = 0;
            dt = 0.002;
            timeSim = 0.0;
            robotIdx = 1;
            %[q,racket,robotIdx] = obj.play(r1,robotIdx,tennisBall.pos);
            
            while timeSim < 4.0 %max(points1,points) < 11 || abs(points1-points2) < 2

                %ballObs = obj.getBallObs();
                %racket = obj.play(r1,ballObs);
                [q,racket,robotIdx] = obj.play(r1,robotIdx,tennisBall.pos);
                %[points1,points2] = obj.judge(points1,points2);
                
                tennisBall.evolve(dt,racket);
                obj.updateAnimation(q);
                timeSim = timeSim + dt;
            end
            
            hold off;
        end
        
        %% FILTERING FOR ROBOTS
        
        function initFilter(obj)
            
            % Initialize EKF
            dim = 6;
            eps = 1e-6; %1e-3;
            C = [eye(3),zeros(3)];
            tennisBall = obj.ball;

            params.C = tennisBall.C;
            params.g = tennisBall.g;
            params.zTable = tennisBall.TABLE.Z;
            params.yNet = tennisBall.NET.Y;
            params.table_length = tennisBall.TABLE.LENGTH; 
            params.table_width = tennisBall.TABLE.WIDTH;
            % coeff of restitution-friction vector
            params.CFTX = tennisBall.TABLE.K(1); 
            params.CFTY = tennisBall.TABLE.K(2); 
            params.CRT = -tennisBall.TABLE.K(3);
            params.ALG = 'RK4'; %'Euler'

            ballFlightFnc = @(x,u,dt) discreteBallFlightModel(x,dt,params);
            % very small but nonzero value for numerical stability
            mats.O = eps * eye(dim);
            mats.C = C;
            mats.M = eps * eye(3);
            obj.filter1 = EKF(dim,ballFlightFnc,mats);
            obj.filter1.initState([tennisBall.pos(:); tennisBall.vel(:)],eps);
            
        end
        
        %% ROBOT STRATEGIES FOR PLAYING TABLE TENNIS 
        
        % FOR NOW WORKING ON ONE ROBOT ONLY        
        function [q,racketStr,robotIdx] = play(obj,robot,robotIdx,ballObs)
            
            % Finite State Automaton
            WAIT = 0;
            PREDICT = 1;
            HIT = 2;            
            stage = obj.FSA;
            
            %loadTennisTableValues();
            dt = 0.002;
            q = obj.q0;
            qs(:,1) = obj.q0;
            qds(:,1) = zeros(7,1);
            [x,xd,o] = robot.calcRacketState([qs;qds]);
            x(:,1) = x; xd(:,1) = xd; o(:,1) = o;
            % FOR DEBUGGING
            ball = obj.ball;
            % land the ball on the centre of opponents court
            desBall(1) = 0.0;
            desBall(2) = obj.NET.Y - obj.TABLE.LENGTH/4;
            desBall(3) = obj.TABLE.Z;
            time2reach = 0.5; % time to reach desired point on opponents court
            
            % if movement finished revert to waiting
            if robotIdx > size(qs,2)
                robotIdx = 1;
                stage = WAIT;
            end
            
            % Estimate the ball state
            filter = obj.filter1;
            filter.linearize(dt,0);
            filter.predict(dt,0);
            filter.update(ballObs,0);            

            time2PassTable = 1.0;
            maxTime2Hit = 0.6;
            maxPredictHorizon = 0.8;
            
            % If it is coming towards the robot consider moving
            velEst = filter.x(4:6);
            if velEst(2) > 0 && stage == WAIT
                stage = PREDICT;
            end
            
            % Estimate ball state and predict ball reaching time
            if stage == PREDICT
                xSave = filter.x;
                PSave = filter.P;
                % update the time it takes to pass table
                yBallEst = filter.x(2);
                tPredIncrement = dt;
                time2PassTable = 0.0;
                while yBallEst <= obj.NET.Y + obj.TABLE.LENGTH/2
                     %filter.linearize(tPredIncrement,0);
                     filter.predict(tPredIncrement,0);
                     yBallEst = filter.x(2);
                     time2PassTable = time2PassTable + tPredIncrement;
                end
                % revert back to saved state
                filter.initState(xSave,PSave);
            end
            
            % Generate trajectory and move
            if stage == PREDICT && time2PassTable <= maxTime2Hit
                
                stage = HIT;
                predictHorizon = maxPredictHorizon;
                predictLen = floor(predictHorizon / dt);
                ballPred = zeros(6,predictLen);
                % FOR DEBUGGING
                filter.initState([ball.pos;ball.vel],PSave);
                for j = 1:predictLen
                    %filter.linearize(dt,0);
                    filter.predict(dt,0);
                    ballPred(:,j) = filter.x;
                end

                % for now only considering the ball positions after table
                ballTime = (1:predictLen) * dt; 

                % Calculate ball outgoing velocities attached to each ball pos
                %%{
                tic
                fast = true; % compute outgoing vel with linear model for speed
                for j = 1:size(ballPred,2)

                    velOut(:,j) = calcBallVelOut3D(desBall,ballPred(1:3,j),time2reach,fast);              
                    % Use the inverse contact model to compute racket vels and normal
                    % at every point                
                    [rp,rv,ro] = calcDesRacketState(ballPred(1:3,j),velOut(:,j),ballPred(4:6,j));
                    racketDes.time(j) = ballTime(j);
                    racketDes.pos(:,j) = rp;
                    racketDes.normal(:,j) = ro;
                    racketDes.vel(:,j) = rv;

                end
                elapsedTimeForCalcDesRacket = toc;
                fprintf('Elapsed time for racket computation: %f sec.\n',...
                    elapsedTimeForCalcDesRacket);
                %}

                % COMPUTE TRAJECTORY HERE

                % define virtual hitting plane (VHP)
                %VHP = -0.6;
                %[qs,qds,qdds] = obj.generate3DVHPTraj(robot,VHP,ballPred,ballTime,q0);
                [qs,qds,qdds] = obj.genOptimalTraj(robot,racketDes,ballPred,ballTime,q0);            
                [qs,qds,qdds] = robot.checkJointLimits(qs,qds,qdds);
                [x,xd,o] = robot.calcRacketState([qs;qds]);

                % Debugging the trajectory generation 
                %h3 = scatter3(x(1,:),x(2,:),x(3,:));
                %h2 = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:));      

            end
            
            % Move the robot
            if stage == HIT 
                q = qs(:,robotIdx);
                racketStr.pos = x(:,robotIdx);
                racketStr.vel = xd(:,robotIdx);
                racketRot = quat2Rot(o(:,robotIdx));
                racketStr.normal = racketRot(:,3);
                robotIdx = robotIdx+1;
            else
                robotIdx = 1;
            end               
            
            obj.FSA = stage;
        end
        
        %% TRAJECTORY GENERATION FOR TABLE TENNIS
        
        % Generate optimal tt trajectories
        function [q,qd,qdd] = genOptimalTraj(obj,robot,racket,ballPred,ballTime,q0)
                  
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
            [qf,qfdot,time2hit] = calcOptimalPoly(robot,racket,ballTime,ballPred,q0);
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
        function [q,qd,qdd] = generate3DVHPTraj(obj,robot,VHP,ballPred,ballTime,q0)
                        
            loadTennisTableValues();
            dof = length(q0);
            time2reach = 0.5; % time to reach desired point on opponents court
            % land the ball on the centre of opponents court
            ballDes(1) = 0.0;
            ballDes(2) = dist_to_table - 3*table_y/2;
            ballDes(3) = table_z + ball_radius;
            [qf,qfdot,timeAtVHP] = calcPolyAtVHP(robot,VHP,time2reach,ballDes,ballPred,ballTime,q0);

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
        function [q,qd,qdd] = generate2DTTTwithVHP(obj,robot,ballPred,ballTime,q0)
                        
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
                qf = robot.invKinematics(R'*racketPos,phiVHP);
                robot.calcJacobian(qf);
                qfdot = robot.jac \ (R'*racketVel);
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
        
        %% Animation functions here
        function initAnimation(obj,q0)           

            % Prepare the animation
            loadTennisTableValues();
            wam = obj.robot1;
            b = obj.ball;

            % get joints, endeffector pos and orientation
            [joint,ee,racket] = wam.drawPosture(q0);

            figure;
            %uisetcolor is useful here
            orange = [0.9100 0.4100 0.1700];
            gray = [0.5020    0.5020    0.5020];
            ballSurfX = b.pos(1) + b.MESH.X;
            ballSurfY = b.pos(2) + b.MESH.Y;
            ballSurfZ = b.pos(3) + b.MESH.Z;
            % transform into base coord.
            h.ball = surf(ballSurfX,ballSurfY,ballSurfZ);
            set(h.ball,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
            hold on;
            h.robot.joints = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
            endeff = [joint(end,:); ee];
            h.robot.endeff = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
            h.robot.racket = fill3(racket(1,:), racket(2,:), racket(3,:), 'r');

            obj.handle = h;
            
            title('Ball-robot interaction');
            grid on;
            axis equal;
            xlabel('x');
            ylabel('y');
            zlabel('z');
            tol_x = 0.1; tol_y = 0.1; tol_z = 0.3;
            xlim([-table_x - tol_x, table_x + tol_x]);
            ylim([dist_to_table - table_length - tol_y, tol_y]);
            zlim([table_z - tol_z, table_z + 2*tol_z]);
            legend('ball','robot');
            fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
            fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
            
        end
        
        % update the animation for robot and the ball
        function updateAnimation(obj,q)
            
            % ANIMATE BOTH THE ROBOT AND THE BALL
            b = obj.ball;
            rob = obj.robot1;
            hball = obj.handle.ball;
            hrobotJ = obj.handle.robot.joints;
            hrobotE = obj.handle.robot.endeff;
            hrobotR = obj.handle.robot.racket;
            [joint,ee,racket] = rob.drawPosture(q);
            endeff = [joint(end,:);ee];
    
            % ball
            set(hball,'XData',b.pos(1) + b.MESH.X);
            set(hball,'YData',b.pos(2) + b.MESH.Y);
            set(hball,'ZData',b.pos(3) + b.MESH.Z);

            % robot joints
            set(hrobotJ,'XData',joint(:,1));
            set(hrobotJ,'YData',joint(:,2));
            set(hrobotJ,'ZData',joint(:,3));
            % robot endeffector
            set(hrobotE,'XData',endeff(:,1));
            set(hrobotE,'YData',endeff(:,2));
            set(hrobotE,'ZData',endeff(:,3));
            % robot racket
            set(hrobotR,'XData',racket(1,:));
            set(hrobotR,'YData',racket(2,:));
            set(hrobotR,'ZData',racket(3,:));

            drawnow;
            %pause(0.001);
            
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end