%% Table Tennis ball class for simulating all sorts of ball interactions

classdef Ball3D < handle
   
    % Fields necessary for all ball models
    properties        
        
        % initial distribution of the ball
        distr
        % air drag
        C
        % gravity
        g
        % position
        pos
        % velocity
        vel
        % radius of the ball
        radius
        % integration method, i.e. RK4 or Euler
        int
        % table structure containing table tennis-related values
        TABLE
        % floor structure containing only z value
        FLOOR
        % net structure
        NET
        % racket structure
        RACKET
        % for drawing
        MESH
        % to give landing information
        isLANDED
    end
    
    % methods that all balls share
    methods (Access = public)   
        
        %% Initialize ball on the ball gun
        % Randomness is due to initial standard deviation initStd
        function obj = Ball3D(initial_distribution)
            
            loadTennisTableValues();   

            obj.C = Cdrag;
            obj.g = gravity;
            obj.radius = ball_radius;
            obj.int = 'RK4';
            
            % table related values
            obj.TABLE.Z = table_z + ball_radius;
            obj.TABLE.WIDTH = table_width;
            obj.TABLE.LENGTH = table_length;
            obj.TABLE.DIST2ROBOT1 = dist_to_table;
            % coeff of restitution-friction vector
            obj.TABLE.K = [CFTX; CFTY; -CRT];
            
            % floor
            obj.FLOOR.Z = floor_level;
            
            % net y value and coeff of restitution of net
            obj.NET.Y = dist_to_table - table_y;
            obj.NET.Xmax = table_width/2 + net_overhang;
            obj.NET.Zmax = table_z + net_height;
            obj.NET.CRN = net_restitution;
            
            % racket restitution
            obj.RACKET.CRR = CRR;
            obj.RACKET.R = racket_radius;
            
            % initialize ball state
            obj.distr = initial_distribution;
            obj.resetState();
            
            % mesh for drawing the ball
            numPoints = 100;
            [ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
            ballMeshX = ball_radius * ballMeshX;
            ballMeshY = ball_radius * ballMeshY;
            ballMeshZ = ball_radius * ballMeshZ;
            obj.MESH.X = ballMeshX;
            obj.MESH.Y = ballMeshY;
            obj.MESH.Z = ballMeshZ;

        end
        
        % Initialize the ball according to the given distribution
        function resetState(obj)
                        
            switch obj.distr.type
                case 'normal'
                    std = obj.distr.init.std;
                    b0 = initBallGun(std,randi(3));
                case 'empirical'
                    b0 = drawSimilarSample(obj.distr.data);
                case 'landing' % draw according to des landing distr
                    b0 = obj.sampleFromLandingDistr();
                case 'workspace' % draw acc. to des workspace distr
                    b0 = obj.sampleFromWorkspaceDistr();
                otherwise
                    error('Distribution not supported');
            end
                
            obj.pos = b0(1:3);
            obj.vel = b0(4:6);
            obj.isLANDED = false;
        end
        
        % sample the ball based on a desired landing distribution
        % sample landing time 
        % sample initial ball to find feasible ball velocities
        % such that ball is above net when it passes
        function b0 = sampleFromLandingDistr(obj)
            
            mean_land_time = 0.5;
            s2_land_time = 0.0004;
            
            par.fast = false;
            par.g = obj.g;
            par.Cdrag = obj.C;
            
            mean_init_pos = obj.distr.init.mean;
            var_init_pos = obj.distr.init.cov;
            mean_land = obj.distr.land.mean;
            mult_land = obj.distr.land.scale;
            znet = -Inf;
            std_init_pos = 0.0;
            std_land = 0.0;
            
            if det(var_init_pos) > 0.0
                std_init_pos = chol(var_init_pos);
            end
           
            
            while znet <= obj.NET.Zmax
                binit = mean_init_pos + std_init_pos * randn(3,1);
                bland = mean_land + mult_land * (rand(2,1)-0.5);
                bland(3) = obj.TABLE.Z + obj.radius;
                tland = mean_land_time + sqrt(s2_land_time) * randn;
                vinit = calcBallVelOut3D(bland,binit,tland,par);
                tnet = abs(obj.NET.Y - binit(2))/vinit(2);
                znet = binit(3) + vinit(3)*tnet + 0.5*obj.g*tnet^2;
            end
                
            b0 = [binit(:);vinit(:)];
        end
        
        % Sample the ball from a desired workspace distr
        % sample time2reach desired workspace point and 
        % sample initial ball to find feasible ball velocities
        % such that ball is above net when it passes and lands on the
        % robots court
        function b0 = sampleFromWorkspaceDistr(obj)
            
            par.fast = false;
            par.g = obj.g;
            par.Cdrag = obj.C;
            par.zTable = obj.TABLE.Z;
            par.radius = obj.radius;
            par.bounce = true;
            par.coeff_bounce = obj.TABLE.K;
            par.dist_to_table = obj.TABLE.DIST2ROBOT1;
            
            mean_reach_time = 1.0;
            s2_reach_time = 0.01;
            mean_init_pos = obj.distr.init.mean;
            var_init_pos = obj.distr.init.cov;
            mean_end_pos = obj.distr.workspace.mean;
            var_end_pos = obj.distr.workspace.cov;
            znet = -Inf;
            xland = Inf;
            yland = Inf;
            y_robot_table_center = obj.TABLE.DIST2ROBOT1 - obj.TABLE.LENGTH/4;
            zland = obj.TABLE.Z + obj.radius;
            
            while znet < obj.NET.Zmax || ...
                  abs(xland) > obj.TABLE.WIDTH/2 || ...
                  abs(yland - y_robot_table_center) > obj.TABLE.LENGTH/4
                 binit = mean_init_pos + chol(var_init_pos) * randn(3,1);
                 bend = mean_end_pos + chol(var_end_pos) * randn(3,1);
                 tend = mean_reach_time + sqrt(s2_reach_time) * randn;
                 vinit = calcBallVel2ReachWorkspace(ballDes,ballPos,time2reach,par);
                 tnet = abs(obj.NET.Y - binit(2))/vinit(2);
                 znet = binit(3) + vinit(3)*tnet - 0.5*obj.g*tnet^2;
                 a = 1/2*obj.g; b = vinit(3); c = binit(3) - zland;
                 tland1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
                 tland2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
                 % pick the one that is positive
                 tland = max(tland1,tland2);
                 assert(tland > 0, 'landing time must be positive!');
                 xland = binit(1) + vinit(1)*tland;
                 yland = binit(2) + vinit(2)*tland;                 
            end
            b0 = [binit(:);vinit(:)];
        end
            
        
        %% Flight models  
        
        % 3D flight model
        function xddot = ballFlightModel(obj,xdot)

            v = sqrt(xdot(1)^2 + xdot(2)^2 + xdot(3)^2);
            xddot(1) = -obj.C * v * xdot(1);
            xddot(2) = -obj.C * v * xdot(2);
            xddot(3) = obj.g - obj.C * v * xdot(3);

            xddot = xddot(:);
        end
        
        %% Integrate the ball state
        
        % Symplectic Integration 
        function evolve(obj,dt,racket)            
            
            % integrate ball flight
            if strcmp(obj.int,'Euler')
                acc = obj.ballFlightModel(obj.vel);
                velNext = obj.vel + dt * acc;
                posNext = obj.pos + dt * velNext;
            elseif strcmp(obj.int,'RK4')
                x = [obj.pos; obj.vel];
                ballFlightFnc = @(x) [x(4:6);obj.ballFlightModel(x(4:6))];
                k1 = dt * ballFlightFnc(x);
                x_k1 = x + k1/2;
                k2 = dt * ballFlightFnc(x_k1);
                x_k2 = x + k2/2;
                k3 = dt * ballFlightFnc(x_k2);
                x_k3 = x + k3;
                k4 = dt * ballFlightFnc(x_k3);
                xNext = x + (k1 + 2*k2 + 2*k3 + k4)/6;
                posNext = xNext(1:3);
                velNext = xNext(4:6);
            else
                error('Not implemented!');
            end
            
            % Check for contact            
            % contact models modify the future states if contact occurs
            [posNext,velNext] = obj.checkContact(dt,posNext,velNext,racket);
            obj.vel = velNext;
            obj.pos = posNext;
            
        end
        
        % Check if dt seconds into the future a contact event occurs
        function [pos,vel] = checkContact(obj,dt,pos,vel,racket)
           
            [pos,vel] = obj.checkContactTable(dt,pos,vel);
            [pos,vel] = obj.checkContactRacket(dt,pos,vel,racket);
            [pos,vel] = obj.checkContactNet(dt,pos,vel);
            [pos,vel] = obj.checkContactFloor(dt,pos,vel);            
            
        end
        
        
        %% Contact models
        
        % K is the coefficient values in y-z directions
        function vel = reboundModel(obj,vel)

            M = diag(obj.TABLE.K);
            vel = M * vel;
        end
        
        % Racket contact model
        function velOut = racketContactModel(obj,velIn,racket)
            % Change ball velocity based on contact model
            racketNormal = racket.normal;
            racketVel = racket.vel;
            
            speedInAlongNormal = racketNormal'*velIn;
            speedRacketAlongNormal = racketNormal'*racketVel;
            velInAlongNormal = speedInAlongNormal.*racketNormal;
            velRacketAlongNormal = speedRacketAlongNormal.*racketNormal;
            % this is kept the same in mirror law
            velInAlongRacket = velIn - velInAlongNormal;
            velOutAlongNormal = velRacketAlongNormal + ...
                obj.RACKET.CRR * (velRacketAlongNormal - velInAlongNormal);
            velOut = velOutAlongNormal + velInAlongRacket;
        end
        
        % rebound for both on and under the table
        % check if rebound event happens within dt seconds
        function [posNext,velNext] = checkContactTable(obj,dt,posNext,velNext)

            % if ball would cross the table
            % consider rebound only            
            cross = sign(posNext(3) - obj.TABLE.Z) ~= ...
                         sign(obj.pos(3) - obj.TABLE.Z);
                     
            % condition for bouncing      
            if cross && abs(posNext(2) - obj.NET.Y) < obj.TABLE.LENGTH/2
                tol = 1e-4;
                dt1 = 0;
                dt2 = dt;
                xBounce = [obj.pos; obj.vel];
                dtBounce = 0.0;
                % doing bisection to find the bounce time
                while abs(xBounce(3) - obj.TABLE.Z) > tol
                    dtBounce = (dt1 + dt2) / 2;
                    xBounce(4:6) = obj.vel + dtBounce * obj.ballFlightModel(obj.vel);
                    xBounce(1:3) = obj.pos + dtBounce * xBounce(4:6);
                    if sign(posNext(3) - obj.TABLE.Z)*(xBounce(3) - obj.TABLE.Z) < 0
                        % increase the time
                        dt1 = dtBounce;
                    else
                        dt2 = dtBounce;
                    end
                end
                % did the ball land on opponents court?
                if ~obj.isLANDED && xBounce(2) < obj.NET.Y
                    obj.isLANDED = true;
                    fprintf('Ball landed at x = %f, y = %f.\n', xBounce(1),xBounce(2));
                end
                % rebound
                xBounce(4:6) = obj.reboundModel(xBounce(4:6));
                % integrate for the rest
                dt = dt - dtBounce;
                velNext = xBounce(4:6) + dt * obj.ballFlightModel(xBounce(4:6));
                posNext = xBounce(1:3) + dt * velNext;
            end
        end
        
        % check contact with racket
        % assumes racket normal stays constant throughout
        % TODO: sometimes does not terminate!
        function [posNext,velNext] = checkContactRacket(obj,dt,posNext,velNext,racket)

            racketPos = racket.pos;
            racketVel = racket.vel;
            racketNormal = racket.normal;
            racketRadius = obj.RACKET.R;
            tol = obj.radius;

            diff = obj.pos - racketPos;
            distToRacketPlane = racketNormal'*diff;
            distOnRacketPlane = sqrt(diff'*diff - distToRacketPlane^2);
            
            % consider contact only  
            % if ball would cross the racket 
            racketPosNext = racketPos + dt * racketVel;
            diffNext = posNext - racketPosNext;
            %cross = sign(diffNext'*racketNormal) ~= ...
            %             sign(diff'*racketNormal);
                     
            distNextToRacketPlane = racketNormal'*diffNext;

            % find the precise hitting time with bisection  
            if distOnRacketPlane < racketRadius && distNextToRacketPlane < obj.radius %cross
                
                tol = 1e-4; % obj.radius
                dt1 = 0;
                dt2 = dt;
                xContact = [obj.pos; obj.vel];
                dtContact = 0.0;
                iter = 0;
                % doing bisection to find the bounce time
                while iter < 5 %abs(distToRacketPlane) > tol
                    dtContact = (dt1 + dt2) / 2;
                    xContact(4:6) = obj.vel + dtContact * obj.ballFlightModel(obj.vel);
                    xContact(1:3) = obj.pos + dtContact * xContact(4:6);
                    racketContactPos = racketPos + dtContact * racketVel;
                    diff = xContact(1:3) - racketContactPos;
                    distToRacketPlane = racketNormal'*diff;
                    if distToRacketPlane > 0
                        % increase the time
                        dt1 = dtContact;
                    else
                        dt2 = dtContact;
                    end
                    iter = iter + 1;
                end
                % racket contact
                xContact(4:6) = obj.racketContactModel(xContact(4:6),racket);
                % integrate for the rest
                dt = dt - dtContact;
                velNext = xContact(4:6) + dt * obj.ballFlightModel(xContact(4:6));
                posNext = xContact(1:3) + dt * velNext;                
                              
                fprintf('Hit at x = %f, y = %f, z = %f\n',...
                     xContact(1),xContact(2),xContact(3));
            end
        end
        
        % If on the other side of the net after integration
        % apply super simplistic net contact model
        function [posNext,velNext] = checkContactNet(obj,dt,posNext,velNext)
            
            tol = obj.radius;
            z = posNext(3);
            y = posNext(2);
            xabs = abs(posNext(1));
            if abs(y - obj.NET.Y) <= tol && z <= obj.NET.Zmax && xabs <= obj.NET.Xmax
                velNext(2) = -obj.NET.CRN * obj.vel(2);                    
            end
            
        end       
        
        % Checks contact, zeros the velocities
        function [pos,vel] = checkContactFloor(obj,dt,pos,vel)
            
            if pos(3) <= obj.FLOOR.Z
                pos(3) = obj.FLOOR.Z;
                vel = zeros(3,1);
            end
        end
        
    end
    
end