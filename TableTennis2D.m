%% 2D Table Tennis for studying trajectory generation

classdef TableTennis2D < handle
    
    properties
        
        % planning related parameters are stored here
        plan
        % table related parameters
        table
        % net is useful for strategies
        net
        % ball class
        ball
        % robot 1
        robot
        % handle structure for drawing animation
        handle
        % draw flag (structure)
        draw
        % noise models (as a structure)
        noise
        % build offline policy (a structure)
        offline
        
    end
    
    methods
        
        %% CONSTRUCTOR AND SUBMETHODS
        function obj = TableTennis2D(rob,q0,opt)
            
            % initialize the robot
            obj.robot = rob;            
            % initialize camera noise
            obj.noise.camera.cov = opt.camera.cov;            
            % initialize a ball    
            obj.ball = Ball2D(opt.distr); 
            % choose method to use for generating trajectories                        
            obj.plan.vhp.flag = opt.plan.vhp.flag;
            obj.plan.vhp.y = opt.plan.vhp.y;
            
            obj.reset_plan(q0);
            obj.init_table();
            obj.init_lookup(opt);
            obj.init_handle(opt,q0);   
            
        end    
        
        % reset planning
        function reset_plan(obj,q0)
            WAIT = 0;
            obj.plan.stage = WAIT;
            obj.plan.idx = 1;
            obj.plan.q = q0;
            obj.plan.qd = zeros(length(q0),1);
        end        

        % initialize table parameters
        function init_table(obj)
            
            loadTennisTableValues();
            % table related values
            obj.table.Z = table_z;
            obj.table.LENGTH = table_length;
            obj.table.DIST = dist_to_table;
            % coeff of restitution-friction vector
            obj.table.K = [CFTY; -CRT];
            
            % net y value and coeff of restitution of net
            obj.net.Y = dist_to_table - table_y;
            obj.net.Zmax = table_z + net_height;
            obj.net.CRN = net_restitution;
        end    
        
        % init handle for recording video and drawing
        function init_handle(obj,opt,q0)            
            % initialize animation
            obj.draw.flag = opt.draw;
            obj.draw.rotate = opt.rotate;
            if opt.draw                
                obj.initAnimation(q0);
                obj.handle.record = false;
                if opt.record
                    obj.handle.record = true;
                    filename = sprintf('tableTennis2DSim%d.avi',randi(100));
                    obj.handle.recordFile = VideoWriter(filename);
                    open(obj.handle.recordFile);
                end
            else
                obj.handle = [];
                obj.handle.record = false;
            end
        end        
        
        % initialize the lookup table parameters
        function init_lookup(obj,opt)
            
            % shall we train an offline lookup table
            obj.offline.train = opt.train;
            obj.offline.use = opt.lookup.flag;
            obj.offline.mode = opt.lookup.mode;
            obj.offline.savefile = opt.lookup.savefile;
            obj.offline.X = [];
            obj.offline.Y = [];
            
            if obj.offline.use || obj.offline.train
                try
                    % load the savefile
                    load(obj.offline.savefile,'X','Y');
                    obj.offline.X = X;
                    obj.offline.Y = Y;
                    obj.offline.B = X \ Y;
                    %obj.train_gp(X,Y);                    
                catch
                    warning('No lookup table found!');
                    obj.offline.X = [];
                    obj.offline.Y = [];
                end
            end            
        end        

        %% MAIN LOOP
        
        % first robot practices 
        function numLands = practice(obj,q0,numTimes)
        
            numLands = 0;
            eps = obj.noise.camera.cov;
            maxSimTime = 3.0;
            dt = 0.01;
            
            % Init filter
            filter = obj.initFilter();
            
            for i = 1:numTimes                                               
                % reset the ball state
                obj.ball.resetState();
                % initialize filter state
                filter.initState([obj.ball.pos;obj.ball.vel],eps);
                % play one turn
                obj.play(dt,q0,filter,maxSimTime);                
                % check landing
                if obj.ball.isLANDED
                    numLands = numLands + 1;
                    % build offline policy
                    if obj.offline.train                        
                        obj.offline.X = [obj.offline.X; obj.offline.b0];
                        obj.offline.Y = [obj.offline.Y; obj.offline.xf];
                    end                    
                end
                fprintf('Iteration: %d\n', i);
                obj.reset_plan(q0);
            end
            
            fprintf('Landed %d/%d.\n',numLands,numTimes);
            if obj.offline.train && ~obj.offline.use
                % save them
                X = obj.offline.X;
                Y = obj.offline.Y;
                save(obj.offline.savefile,'X','Y');
            end
            
            % make sure recording is closed
            if obj.handle.record
                close(obj.handle.recordFile);
            end
        end   
        
        % first robot plays solo once
        function play(obj,dt,q0,filter,timeMax)
            
            timeSim = 0.0;
            % initialize q and x
            qd0 = zeros(length(q0),1);
            rotAngle = obj.draw.rotate;
            [x,xd,o] = obj.robot.getEndEffectorState(q0,qd0,rotAngle);

            while timeSim < timeMax      
                
                % evolve ball according to racket and get estimate
                filter = obj.getBallEstimate(dt,filter,x,xd,o);
                [q,qd] = obj.planFiniteStateMachine(filter,q0,dt);
                [x,xd,o] = obj.robot.getEndEffectorState(q,qd,rotAngle);
                
                if obj.draw.flag
                    obj.updateAnimation(q);
                end
                timeSim = timeSim + dt;                
            end
            
            % clear the ball path predicted and robot generated traj
            if obj.draw.flag && isfield(obj.handle,'ballPred')
                set(obj.handle.robotCartesian,'Visible','off');
                set(obj.handle.ballPred,'Visible','off');
            end
        end        
        
        % planning using a Virtual planning plane (VPP)
        % generally over the net
        % using a Finite State Machine to plan when to hit/stop
        function [q,qd] = planFiniteStateMachine(obj,filter,q0,dt)
            
            WAIT = 0;
            PREDICT = 1;
            HIT = 2;
            FINISH = 3; % only when practicing solo
            
            % If it is coming towards the robot consider predicting
            velEst = filter.x(3:4);
            if velEst(1) > 0 && obj.plan.stage == WAIT
                obj.plan.stage = PREDICT;
            end

            table_center = obj.table.DIST - obj.table.LENGTH/2;
            %if stage == PREDICT && time2Passtable <= minTime2Hit       
            if obj.plan.stage == PREDICT && filter.x(1) > table_center && filter.x(3) > 0.5   
                predictTime = 1.2;
                [ballPred,ballTime,numBounce,time2Passtable] = ...
                    predictBallPath2D(dt,predictTime,filter,obj.table);
                if checkBounceOnOppTable2D(filter,obj.table)
                    obj.plan.stage = FINISH;
                elseif numBounce ~= 1
                    disp('Ball does not bounce once! Not hitting!');
                    obj.plan.stage = FINISH;
%                 elseif ~checkIfBallIsInsideWorkspace(obj.robot,ballPred)
%                     disp('No intersection with workspace! Not hitting!');
%                     obj.plan.stage = FINISH;
                else
                    obj.plan.idx = 0;
                    obj.plan.stage = HIT;
                    % If we're training an offline model save optimization result
                    if obj.offline.train || obj.offline.use
                        obj.offline.b0 = filter.x';
                    end
                    obj.returnBall2Center(ballPred,ballTime,q0,dt);            
                end

            end % end predict      

            % Move the robot
            if obj.plan.stage == HIT 
                obj.plan.idx = obj.plan.idx+1;
            end
            % if movement finished revert to waiting
            if obj.plan.idx > size(obj.plan.q,2)
                obj.plan.idx = size(obj.plan.q,2);
                obj.plan.stage = FINISH;
            end
            
            q = obj.plan.q(:,obj.plan.idx);
            qd = obj.plan.qd(:,obj.plan.idx);
        end  
        
        % Fix a desired landing point and desired landing time
        % and calculate racket variables over ball estimated trajectory
        function racketDes = planRacket(obj,ballDes,ballPred,ballTime,time2reach,q0)
            
            %Calculate ball outgoing velocities attached to each ball pos
            fast = false;
            racketDes = calcRacketStrategy2D(ballDes,ballPred,ballTime,time2reach,fast);
            
            % Initialize solution for optimal poly
            timeEst = 0.8;
            q0dot = zeros(length(q0),1);
            x0 = [q0;q0dot;timeEst];
            racketDes.est = x0;
        end  
        
        % Loads lookup table parameters by finding the 
        % closest table ball-estimate entry 
        function [qf,qfdot,T] = lookup(obj)            
       
            switch obj.offline.mode
                case 'regress'
                    val = obj.offline.b0 * obj.offline.B;
                case 'closest'      
                    Xs = obj.offline.X;
                    Ys = obj.offline.Y;
                    bstar = obj.offline.b0;
                    N = size(Xs,1);
                    % find the closest point among Xs
                    diff = repmat(bstar,N,1) - Xs;
                    [~,idx] = min(diag(diff*diff'));
                    val = Ys(idx,:);              
                otherwise
                    error('lookup mode not supported!');
            end
            dofs = (length(val) - 1) / 2;
            qf = val(1:dofs)';
            qfdot = val(dofs+1:2*dofs)';
            T = val(end);
        end        
        
        % Fix a desired landing point and desired landing time
        % as well as a desired return time to q0
        % For now only two methods : VHP and free Time
        function returnBall2Center(obj,ballPred,ballTime,q0,dt)                
            
            dofs = length(q0);
            rotAngle = obj.draw.rotate;
            % land the ball on the centre of opponents court
            ballDes(1) = obj.table.DIST - 3*obj.table.LENGTH/4;
            ballDes(2) = obj.table.Z;   
            q0dot = zeros(dofs,1);
            time2return = 1.0; % time for robot to go back to q0 after hit
            
            if obj.offline.use
                [qf,qfdot,T] = obj.lookup();
            else
                time2reach = 0.8; % time to reach desired point on opponents court    

                % Compute traj here
                if obj.plan.vhp.flag
                    [qf,qfdot,T] = calcPolyAtVHP2D(obj.robot,obj.plan.vhp.y,time2reach,ballDes,ballPred,ballTime,q0);
                else
                    racketDes = obj.planRacket(ballDes,ballPred,ballTime,time2reach,q0);
                    [qf,qfdot,T] = calcOptimalPoly(obj.robot,racketDes,q0,time2return);
                end
                % If we're training an offline model save optimization result
                if obj.offline.train
                    obj.offline.xf = [qf',qfdot',T];
                end
            end
            
            [q,qd,qdd] = generateSpline(dt,q0,q0dot,qf,qfdot,T,time2return);
            [q,qd,qdd] = obj.robot.checkJointLimits(q,qd,qdd);
            [x,xd,o] = obj.robot.getEndEffectorState(q,qd,rotAngle);

            if obj.draw.flag
                % Debugging the trajectory generation                 
                obj.handle.robotCartesian = scatter(x(1,:),x(2,:));
                obj.handle.ballPred = scatter(ballPred(1,:),ballPred(2,:));
            end
            
            obj.plan.q = q;
            obj.plan.qd = qd;
        end        
  
        
        %% FILTERING FOR ROBOTS TO GET BALL OBSERVATION
        
        function filter = initFilter(obj)
            
            % Initialize EKF
            dim = 4;
            eps = 1e-6; %1e-3;
            C = [eye(2),zeros(2)];
            tennisBall = obj.ball;

            params.C = tennisBall.C;
            params.g = tennisBall.g;
            params.radius = tennisBall.radius;
            params.zTable = obj.table.Z;
            params.yNet = obj.net.Y;
            params.table_length = obj.table.LENGTH; 
            % coeff of restitution-friction vector
            params.CFTY = obj.table.K(1); 
            params.CRT = -obj.table.K(2);
            params.ALG = 'RK4'; 

            ballFlightFnc = @(x,u,dt) discreteBallFlightModel2D(x,dt,params);
            % very small but nonzero value for numerical stability
            mats.O = eps * eye(dim);
            mats.C = C;
            mats.M = eps * eye(2);
            filter = EKF(dim,ballFlightFnc,mats);
            filter.initState([tennisBall.pos(:); tennisBall.vel(:)],eps);            
        end    
        
        function obs = emulateCamera(obj)
            
            if det(obj.noise.camera.cov) > 0
                std = chol(obj.noise.camera.cov);
            else
                std = 0.0;
            end
            obs = obj.ball.pos + std * randn(2,1);
        end
        
        function filter = getBallEstimate(obj,dt,filter,x,xd,o)
            racket.pos = x;
            racket.vel = xd;
            racket.normal = o(1:2,1,3); 
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
            rob = obj.robot;
            ball = obj.ball;

            % get joints, endeffector pos and orientation
            [joint,ee,racket] = rob.drawPosture(q,obj.draw.rotate);
            endeff = [joint(end,:); ee];

            % edit: reduced to half screen size
            scrsz = get(groot,'ScreenSize');
            figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
            %uisetcolor is useful here to determine these 3-vectors
            orange = [0.9100 0.4100 0.1700];
            gray = [0.5020 0.5020 0.5020];
            % transform into base coord.
            h.ball = scatter(ball.pos(1),ball.pos(2),50,orange,'filled');
            hold on;
            h.robot.joints = plot(joint(:,1),joint(:,2),'k','LineWidth',10);
            h.robot.endeff = plot(endeff(:,1),endeff(:,2),'Color',gray,'LineWidth',10);
            h.robot.racket = line(racket(:,1),racket(:,2),'color',[1 0 0],'LineWidth',10);
            obj.handle = h;
            
            grid on;
            axis equal;
            xlabel('y');
            ylabel('z');

            % define table for 2d case
            table_z = floor_level + table_height;
            table_y = table_length/2;
            table_width_2d = 0.05;
            net_width_2d = 0.01;

            T1 = [dist_to_table - table_length; 
                table_z - table_width_2d];
            T2 = [dist_to_table - table_length;
                table_z];
            T3 = [dist_to_table;
                table_z];
            T4 = [dist_to_table;
                table_z - table_width_2d];
            T = [T1,T2,T3,T4]';

            net1 = [dist_to_table - table_y - net_width_2d;
                    table_z];
            net2 = [dist_to_table - table_y + net_width_2d;
                    table_z];
            net3 = [dist_to_table - table_y + net_width_2d;
                    table_z + net_height];
            net4 = [dist_to_table - table_y - net_width_2d;
                    table_z + net_height];

            net = [net1,net2,net3,net4]';

            tol_y = 0.1; tol_z = 0.3;
            xlim([dist_to_table - table_length - tol_y, 3*tol_y]);
            ylim([table_z - 3*tol_z, table_z + 4*tol_z]);
            legend('ball','robot');
            fill(T(:,1),T(:,2),[0 0.7 0.3]);
            fill(net(:,1),net(:,2),[0 0 0]);  
        end        
        
        % update the animation for robot and the ball
        function updateAnimation(obj,q)
            
            % ANIMATE BOTH THE ROBOT AND THE BALL
            b = obj.ball;
            rob = obj.robot;
            hball = obj.handle.ball;
            hrobotJ = obj.handle.robot.joints;
            hrobotE = obj.handle.robot.endeff;
            hrobotR = obj.handle.robot.racket;
            rotangle = obj.draw.rotate;
            [joint,ee,racket] = rob.drawPosture(q,rotangle);
            endeff = [joint(end,:); ee];
    
            % ball
            set(hball,'XData',b.pos(1));
            set(hball,'YData',b.pos(2));

            % robot joints
            set(hrobotJ,'XData',joint(:,1));
            set(hrobotJ,'YData',joint(:,2));
            % robot endeffector
            set(hrobotE,'XData',endeff(:,1));
            set(hrobotE,'YData',endeff(:,2));
            % robot racket
            set(hrobotR,'XData',racket(:,1));
            set(hrobotR,'YData',racket(:,2));

            drawnow;
            %pause(0.001);
            
            if obj.handle.record
                frame = getframe(gcf);
                writeVideo(obj.handle.recordFile,frame);
            end            
        end    
        
    end
end