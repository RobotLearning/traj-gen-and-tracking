%% Table Tennis class for simulating a solo game

classdef TableTennis < handle
    
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
        function obj = TableTennis(wam,q0,opt)
            
            % initialize the robot
            obj.robot = wam;            
            % initialize camera noise
            obj.noise.camera.cov = opt.camera.cov;            
            % initialize a ball    
            obj.ball = Ball(opt.distr); 
            % choose method to use for generating trajectories                        
            obj.plan.vhp.flag = opt.plan.vhp.flag;
            obj.plan.vhp.y = opt.plan.vhp.y;
            
            obj.reset_plan(q0);
            obj.init_lookup(opt);
            obj.init_table(); 
            obj.init_handle(opt,q0);   
            
        end
        
        % reset planning
        function reset_plan(obj,q0)
            WAIT = 0;
            obj.plan.stage = WAIT;
            obj.plan.idx = 1;
            obj.plan.q = q0;
            obj.plan.qd = zeros(7,1);
        end
        
        % initialize table parameters
        function init_table(obj)
            
            loadTennisTableValues();
            % table related values
            obj.table.Z = table_z;
            obj.table.WIDTH = table_width;
            obj.table.LENGTH = table_length;
            obj.table.DIST = dist_to_table;
            % coeff of restitution-friction vector
            obj.table.K = [CFTX; CFTY; -CRT];
            
            % net y value and coeff of restitution of net
            obj.net.Y = dist_to_table - table_y;
            obj.net.Xmax = table_width/2 + net_overhang;
            obj.net.Zmax = table_z + net_height;
            obj.net.CRN = net_restitution;    
        end
        
        % init handle for recording video and drawing
        function init_handle(obj,opt,q0)            
            % initialize animation
            obj.draw.flag = opt.draw;
            if opt.draw                
                obj.initAnimation(q0);
                obj.handle.record = false;
                if opt.record
                    obj.handle.record = true;
                    filename = sprintf('tableTennisSim%d.avi',randi(100));
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
        
        % Train independent GPS
        function train_gp(obj,X,Y)            
            hp.type = 'squared exponential iso';
            hp.l = 1/4;
            hp.scale = 1;
            hp.noise.var = 0.0;
            ndofs = 7;
            num_dims = 2*ndofs + 1;
            for i = 1:num_dims
                gp{i} = GP(hp,X',Y(:,i));
            end
            obj.offline.GP = gp;            
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
            qd0 = zeros(7,1);
            [x,xd,o] = obj.robot.calcRacketState([q0;qd0]);

            while timeSim < timeMax      
                
                % evolve ball according to racket and get estimate
                filter = obj.getBallEstimate(dt,filter,x,xd,o);
                [q,qd] = obj.planFiniteStateMachine(filter,q0,dt);
                [x,xd,o] = obj.robot.calcRacketState([q;qd]);
                
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
        
        %% STRATEGIES/TRAJ GEN FOR TABLE TENNIS
        
        % planning using a Virtual planning plane (VPP)
        % generally over the net
        % using a Finite State Machine to plan when to hit/stop
        function [q,qd] = planFiniteStateMachine(obj,filter,q0,dt)
            
            WAIT = 0;
            PREDICT = 1;
            HIT = 2;
            FINISH = 3; % only when practicing solo
            
            % If it is coming towards the robot consider predicting
            velEst = filter.x(4:6);
            if velEst(2) > 0 && obj.plan.stage == WAIT
                obj.plan.stage = PREDICT;
            end

            table_center = obj.table.DIST - obj.table.LENGTH/2;
            %if stage == PREDICT && time2Passtable <= minTime2Hit       
            if obj.plan.stage == PREDICT && filter.x(2) > table_center && filter.x(5) > 0.5   
                predictTime = 1.2;
                [ballPred,ballTime,numBounce,time2Passtable] = ...
                    predictBallPath(dt,predictTime,filter,obj.table);
                if checkBounceOnOppTable(filter,obj.table)
                    obj.plan.stage = FINISH;
                elseif numBounce ~= 1
                    disp('Ball does not bounce once! Not hitting!');
                    obj.plan.stage = FINISH;
                elseif ~checkIfBallIsInsideWorkspace(obj.robot,ballPred)
                    disp('No intersection with workspace! Not hitting!');
                    obj.plan.stage = FINISH;
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
            racketDes = calcRacketStrategy(ballDes,ballPred,ballTime,time2reach,fast);
            
            % Initialize solution for optimal poly
            timeEst = 0.8;
            q0dot = zeros(7,1);
            x0 = [q0;q0dot;timeEst];
            racketDes.est = x0;
        end
        
        % Loads lookup table parameters by finding the 
        % closest table ball-estimate entry 
        function [qf,qfdot,T] = lookup(obj)            
       
            switch obj.offline.mode
                case 'regress'
                    val = obj.offline.b0 * obj.offline.B;
                    %{
                    numdims = 2*dofs + 1;
                    val = zeros(1,numdims);
                    for i = 1:numdims
                        val(i) = obj.offline.GP{i}.predict(obj.offline.b0);
                    end
                    %}
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
            
            dofs = 7;
            % land the ball on the centre of opponents court
            ballDes(1) = 0.0;
            ballDes(2) = obj.table.DIST - 3*obj.table.LENGTH/4;
            ballDes(3) = obj.table.Z;   
            q0dot = zeros(dofs,1);
            time2return = 1.0; % time for robot to go back to q0 after hit
            
            if obj.offline.use
                [qf,qfdot,T] = obj.lookup();
            else
                time2reach = 0.8; % time to reach desired point on opponents court    

                % Compute traj here
                if obj.plan.vhp.flag
                    [qf,qfdot,T] = calcPolyAtVHP(obj.robot,obj.plan.vhp.y,time2reach,ballDes,ballPred,ballTime,q0);
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
            [x,xd,o] = obj.robot.calcRacketState([q;qd]);

            if obj.draw.flag
                % Debugging the trajectory generation                 
                obj.handle.robotCartesian = scatter3(x(1,:),x(2,:),x(3,:));
                obj.handle.ballPred = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:));
            end
            
            obj.plan.q = q;
            obj.plan.qd = qd;
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
            params.radius = tennisBall.radius;
            params.zTable = obj.table.Z;
            params.yNet = obj.net.Y;
            params.table_length = obj.table.LENGTH; 
            params.table_width = obj.table.WIDTH;
            % coeff of restitution-friction vector
            params.CFTX = obj.table.K(1); 
            params.CFTY = obj.table.K(2); 
            params.CRT = -obj.table.K(3);
            params.ALG = 'RK4'; 

            ballFlightFnc = @(x,u,dt) discreteBallFlightModel(x,dt,params);
            % very small but nonzero value for numerical stability
            mats.O = eps * eye(dim);
            mats.C = C;
            mats.M = eps * eye(3);
            filter = EKF(dim,ballFlightFnc,mats);
            filter.initState([tennisBall.pos(:); tennisBall.vel(:)],eps);            
        end
        
        function obs = emulateCamera(obj)
            
            if det(obj.noise.camera.cov) > 0
                std = chol(obj.noise.camera.cov);
            else
                std = 0.0;
            end
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
            wam = obj.robot;
            ball = obj.ball;

            % get joints, endeffector pos and orientation
            [joint,ee,racket] = wam.drawPosture(q);
            endeff = [joint(end,:); ee];

            % edit: reduced to half screen size
            scrsz = get(groot,'ScreenSize');
            figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
            %uisetcolor is useful here to determine these 3-vectors
            orange = [0.9100 0.4100 0.1700];
            gray = [0.5020 0.5020 0.5020];
            lightgray = [0.8627    0.8627    0.8627];
            white = [0.9412 0.9412 0.9412];
            black2 = [0.3137    0.3137    0.3137];
            red = [1.0000    0.2500    0.2500];
            ballSurfX = ball.pos(1) + ball.MESH.X;
            ballSurfY = ball.pos(2) + ball.MESH.Y;
            ballSurfZ = ball.pos(3) + ball.MESH.Z;
            % transform into base coord.
            h.ball = surf(ballSurfX,ballSurfY,ballSurfZ);
            set(h.ball,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
            hold on;
            h.robot.joints = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
            h.robot.endeff = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
            h.robot.racket = fill3(racket(1,:), racket(2,:), racket(3,:),red);

            obj.handle = h;
            
            %title('Ball-robot interaction');
            grid on;
            axis equal;
            xlabel('x');
            ylabel('y');
            zlabel('z');
            tol_x = 0.2; tol_y = 0.1; tol_z = 0.3;
            xlim([-table_x - tol_x, table_x + tol_x]);
            ylim([dist_to_table - table_length - tol_y, 3*tol_y]);
            zlim([floor_level, table_z + 5*tol_z]);
            %legend('ball','robot');
            %fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
            % faces matrix 6x4
            F = [1 2 3 4;
                 5 6 7 8;
                 1 2 6 5;
                 2 3 7 6;
                 3 4 8 7;
                 4 1 5 8];
            table = patch('Faces',F,'Vertices',T,'FaceColor',[0 0.7 0.3]);
            cdata = [0 0.7 0.3;
                     0 0.7 0.3;
                     repmat(black2,4,1)];
            set(table,'FaceColor','flat','FaceVertexCData',cdata);
            
            % plot the robot stand
            patch('Faces',F,'Vertices',SR,'FaceColor',black2);
            
            % instead draw 14 black thin lines
            numpts = 1000;
            num_horz_lines = 10;
            num_vert_lines = 50;
            tol = 0.02;
            x_nets = repmat(linspace(-table_x-net_overhang,table_x + net_overhang,numpts),num_horz_lines,1);
            y_nets = repmat(linspace(dist_to_table - table_y,dist_to_table - table_y,numpts),num_horz_lines,1);
            z_nets = repmat(linspace(table_z+tol,table_z+net_height-tol,num_horz_lines)',1,numpts);
            plot3(x_nets',y_nets',z_nets','Color',black2,'LineWidth',0.5);
            x_nets = repmat(linspace(-table_x-net_overhang,table_x + net_overhang,num_vert_lines)',1,100);
            y_nets = repmat(linspace(dist_to_table - table_y,dist_to_table - table_y,100),num_vert_lines,1);
            z_nets = repmat(linspace(table_z+tol,table_z+net_height-tol,100),num_vert_lines,1);
            plot3(x_nets',y_nets',z_nets','Color',black2,'LineWidth',0.5);
            topline_x = linspace(-table_x-net_overhang,table_x+net_overhang,numpts);
            topline_y = (dist_to_table-table_y) * ones(1,numpts);
            topline_z = (table_z + net_height) * ones(1,numpts);
            plot3(topline_x,topline_y,topline_z,'Color',white,'LineWidth',2);
            botline_x = linspace(-table_x-net_overhang,table_x+net_overhang,numpts);
            botline_y = (dist_to_table-table_y) * ones(1,numpts);
            botline_z = (table_z+tol) * ones(1,numpts);
            plot3(botline_x,botline_y,botline_z,'k','LineWidth',2);
            lefthang_x = (-table_x-net_overhang) * ones(1,100);
            lefthang_y = (dist_to_table-table_y) * ones(1,100);
            lefthang_z = linspace(table_z+tol,table_z+net_height,100);
            plot3(lefthang_x,lefthang_y,lefthang_z,'k','LineWidth',4);
            righthang_x = (table_x+net_overhang) * ones(1,100);
            plot3(righthang_x,lefthang_y,lefthang_z,'k','LineWidth',4);
            %fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
            
            % three white vert lines and two horizontal lines
            tol = 0.012;
            line1_x = repmat(linspace(-table_x+tol,table_x-tol,3)',1,numpts);
            line1_y = repmat(linspace(dist_to_table-table_length,dist_to_table,numpts),3,1);
            line1_z = repmat(table_z * ones(1,numpts),3,1);
            plot3(line1_x',line1_y',line1_z','w','LineWidth',3);
            line2_x = repmat(linspace(-table_x,table_x,numpts),2,1);
            line2_y = repmat(linspace(dist_to_table-table_length+tol,dist_to_table-tol,2)',1,numpts);
            line2_z = repmat(table_z * ones(1,numpts),2,1);
            plot3(line2_x',line2_y',line2_z','w','LineWidth',5);
            
            % fill also the virtual hitting plane
            if obj.plan.vhp.flag
                V1 = [table_center - table_x; 
                    obj.VHP.Y; 
                    table_z - 2*tol_z];
                V2 = [table_center + table_x;
                    obj.VHP.Y;
                    table_z - 2*tol_z];
                V3 = [table_center + table_x;
                    obj.VHP.Y;
                    table_z + 2*tol_z];
                V4 = [table_center - table_x;
                    obj.VHP.Y;
                    table_z + 2*tol_z];
                V = [V1,V2,V3,V4]';
                text(table_center-table_x,obj.VHP.Y,table_z+2*tol_z-0.05,'VHP')
                fill3(V(:,1),V(:,2),V(:,3),[0 0.7 0.3],'FaceAlpha',0.2);
            end
            
            % change view angle
            az = -81.20; % azimuth
            el = 16.40; % elevation
            % angles manually tuned
            view(az,el);
            
            
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
            
            if obj.handle.record
                frame = getframe(gcf);
                writeVideo(obj.handle.recordFile,frame);
            end

            
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end