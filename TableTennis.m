%% Table Tennis class for simulating a match between two robots!

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
        % handle structure for drawing animation
        handle
        % draw flag
        drawFlag
        % noise models (as a structure)
        noise
        
    end
    
    methods
        
        %% CONSTRUCTOR
        function obj = TableTennis(wam,wam2,q0,std,draw)
            
            % initialize two robots
            obj.robot1 = wam;
            obj.robot2 = wam2;
            
            % initialize the noise models
            obj.noise.ballInit.pos.std = std.pos;
            obj.noise.ballInit.vel.std = std.vel;
            obj.noise.camera.std = std.camera;
            
            % initialize a ball            
            obj.ball = Ball(0.0,0.0); 
                       
            % initialize animation
            obj.drawFlag = draw;
            if draw
                obj.initAnimation(q0);  
            else
                obj.handle = [];
            end
            
            loadTennisTableValues();
            % table related values
            obj.TABLE.Z = table_z;
            obj.TABLE.WIDTH = table_width;
            obj.TABLE.LENGTH = table_length;
            obj.TABLE.DIST = dist_to_table;
            % coeff of restitution-friction vector
            obj.TABLE.K = [CFTX; CFTY; -CRT];
            
            % net y value and coeff of restitution of net
            obj.NET.Y = dist_to_table - table_y;
            obj.NET.Xmax = table_width/2 + net_overhang;
            obj.NET.Zmax = table_z + net_height;
            obj.NET.CRN = net_restitution;           
            
        end
        
        %% MAIN LOOP
        
        % first robot practices 
        function numLands = practice(obj,numTimes)
        
            stdPos = obj.noise.ballInit.pos.std;
            stdVel = obj.noise.ballInit.vel.std;
            numLands = 0;
            eps = 0.0;
            maxSimTime = 3.0;
            
            % Initialize arm posture
            % initialize the arm with zero velocity on the right hand side
            q0 = [1.8; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
            
            % Init filter
            filter = obj.initFilter();
            
            for i = 1:numTimes
                                               
                % initialize a ball
                obj.ball = Ball(stdPos,stdVel);
                % initialize filter state
                filter.initState([obj.ball.pos;obj.ball.vel],eps);
                % play one turn
                obj.play(q0,filter,maxSimTime);
                
                % check landing
                if obj.ball.isLANDED
                    numLands = numLands + 1;
                end
            end
        end
        
        % first robot plays solo once
        function play(obj,q0,filter,timeMax)
            
            r1 = obj.robot1;
            dt = 0.01;
            timeSim = 0.0;
            idx = 1;
            maxTime2Hit = 0.6;
            
            WAIT = 0;
            PREDICT = 1;
            HIT = 2;
            FINISH = 3; % only when practicing solo
            stage = WAIT; 
            
            % initialize q and x
            q = q0; qd0 = zeros(7,1);
            [x,xd,o] = r1.calcRacketState([q0;qd0]);

            while timeSim < timeMax
                
                % evolve ball according to racket and get estimate
                filter = obj.getBallEstimate(dt,filter,x(:,idx),xd(:,idx),o(:,idx));

                % If it is coming towards the robot consider predicting
                velEst = filter.x(4:6);
                if velEst(2) > 0 && stage == WAIT
                    stage = PREDICT;
                end    
                
                if stage == PREDICT
                    [ballPred,ballTime,numBounce,time2PassTable] = predictBallPath(dt,filter,obj.TABLE);
                    if checkBounceOnOppTable(filter,obj.TABLE)
                        stage = FINISH;
                    end
                end
                
                if stage == PREDICT && time2PassTable <= maxTime2Hit         
                                       
                    if numBounce ~= 1
                        disp('Ball is not valid! Not hitting!');
                        stage = FINISH;
                    else
                        idx = 0;
                        stage = HIT;
                        
                        % Calculate ball outgoing velocities attached to each ball pos
                        tic
                        fast = true;
                        % land the ball on the centre of opponents court
                        desBall(1) = 0.0;
                        desBall(2) = obj.TABLE.DIST - 3*obj.TABLE.LENGTH/4;
                        desBall(3) = obj.TABLE.Z;
                        time2reach = 0.5; % time to reach desired point on opponents court
                        racketDes = calcRacketStrategy(desBall,ballPred,ballTime,time2reach,fast);
                        elapsedTimeForCalcDesRacket = toc;
                        fprintf('Elapsed time for racket computation: %f sec.\n',...
                            elapsedTimeForCalcDesRacket);

                        % Compute traj here
                        %VHP = -0.4;
                        %[q,qd,qdd] = generate3DTTTwithVHP(r1,VHP,ballPred,ballTime,q0);
                        [q,qd,qdd] = generateOptimalTTT(r1,racketDes,dt,q0);
                        [q,qd,qdd] = r1.checkJointLimits(q,qd,qdd);
                        [x,xd,o] = r1.calcRacketState([q;qd]);

                        if obj.drawFlag
                            % Debugging the trajectory generation 
                            h3 = scatter3(x(1,:),x(2,:),x(3,:));
                            h2 = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:));
                        end
                    
                    end

                end % end predict      
                
                % Move the robot
                if stage == HIT 
                    idx = idx+1;
                end

                % if movement finished revert to waiting
                if idx > size(q,2)
                    idx = size(q,2);
                    stage = FINISH;
                end
                
                if obj.drawFlag
                    obj.updateAnimation(q(:,idx));
                end
                timeSim = timeSim + dt;
                
            end
            
            if obj.drawFlag && exist('h2')
                set(h2,'Visible','off');
                set(h3,'Visible','off');
            end
        end
        
        %% FILTERING FOR ROBOTS TO GET BALL OBSERVATION
        
        function filter = initFilter(obj)
            
            % Initialize EKF
            dim = 6;
            eps = 1e-6; %1e-3;
            C = [eye(3),zeros(3)];
            tennisBall = obj.ball;

            params.C = tennisBall.C;
            params.g = tennisBall.g;
            params.zTable = obj.TABLE.Z;
            params.yNet = obj.NET.Y;
            params.table_length = obj.TABLE.LENGTH; 
            params.table_width = obj.TABLE.WIDTH;
            % coeff of restitution-friction vector
            params.CFTX = obj.TABLE.K(1); 
            params.CFTY = obj.TABLE.K(2); 
            params.CRT = -obj.TABLE.K(3);
            params.ALG = 'RK4'; %'Euler'

            ballFlightFnc = @(x,u,dt) discreteBallFlightModel(x,dt,params);
            % very small but nonzero value for numerical stability
            mats.O = eps * eye(dim);
            mats.C = C;
            mats.M = eps * eye(3);
            filter = EKF(dim,ballFlightFnc,mats);
            filter.initState([tennisBall.pos(:); tennisBall.vel(:)],eps);
            
        end
        
        function obs = emulateCamera(obj)
            
            std = obj.noise.camera.std;
            obs = obj.ball.pos + std * randn(3,1);
        end
        
        function filter = getBallEstimate(obj,dt,filter,x,xd,o)
            racket.pos = x;
            racket.vel = xd;
            racketRot = quat2Rot(o);
            racket.normal = racketRot(:,3);
            obj.ball.evolve(dt,racket);
            obs = obj.emulateCamera();

            % Estimate the ball state
            filter.linearize(dt,0);
            filter.predict(dt,0);
            filter.update(obs,0);                
        end
        
        %% Animation functions here
        function initAnimation(obj,q)           

            % Prepare the animation
            loadTennisTableValues();
            wam = obj.robot1;
            b = obj.ball;

            % get joints, endeffector pos and orientation
            [joint,ee,racket] = wam.drawPosture(q);
            endeff = [joint(end,:); ee];

            figure;
            %uisetcolor is useful here to determine these 3-vectors
            orange = [0.9100 0.4100 0.1700];
            gray = [0.5020 0.5020 0.5020];
            ballSurfX = b.pos(1) + b.MESH.X;
            ballSurfY = b.pos(2) + b.MESH.Y;
            ballSurfZ = b.pos(3) + b.MESH.Z;
            % transform into base coord.
            h.ball = surf(ballSurfX,ballSurfY,ballSurfZ);
            set(h.ball,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
            hold on;
            h.robot.joints = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
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