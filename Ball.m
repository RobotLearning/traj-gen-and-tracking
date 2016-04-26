%% Table Tennis ball class for simulating all sorts of ball interactions

classdef Ball < handle
   
    % Fields necessary for all ball models
    properties        
        
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
        function obj = Ball(meanInit,varInit)
            
            loadTennisTableValues();
            
            initState = meanInit + chol(varInit) * randn(6,1);
            obj.pos = initState(1:3);
            obj.vel = initState(4:6);
            
            obj.C = Cdrag;
            obj.g = gravity;
            obj.radius = ball_radius;
            obj.int = 'Euler';
            
            % table related values
            obj.TABLE.Z = table_z;
            obj.TABLE.WIDTH = table_width;
            obj.TABLE.LENGTH = table_length;
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
            
            % mesh for drawing the ball
            numPoints = 100;
            [ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
            ballMeshX = ball_radius * ballMeshX;
            ballMeshY = ball_radius * ballMeshY;
            ballMeshZ = ball_radius * ballMeshZ;
            obj.MESH.X = ballMeshX;
            obj.MESH.Y = ballMeshY;
            obj.MESH.Z = ballMeshZ;
            
            obj.isLANDED = false;
        end
        
        %% Flight models  
        % 2D flight model
        function xddot = ballFlightModel2D(obj,xdot)

            v = sqrt(xdot(1)^2 + xdot(2)^2);
            xddot(1) = -obj.C * v * obj.vel(1);
            xddot(2) = obj.g - obj.C * v * obj.vel(2);

            xddot = xddot(:);
        end
        
        % 3D flight model
        function xddot = ballFlightModel3D(obj,xdot)

            v = sqrt(xdot(1)^2 + xdot(2)^2 + xdot(3)^2);
            xddot(1) = -obj.C * v * xdot(1);
            xddot(2) = -obj.C * v * xdot(2);
            xddot(3) = obj.g - obj.C * v * xdot(3);

            xddot = xddot(:);
        end
        
        %% Integrate the ball state
        
        % Symplectic Integration 
        function evolve(obj,dt,racket)            

            acc = obj.ballFlightModel3D(obj.vel);
            %acc = obj.ballFlightModel2D();
            
            % integrate ball flight
            if strcmp(obj.int,'Euler')
                velNext = obj.vel + dt * acc;
                posNext = obj.pos + dt * velNext;
            else
                error('Not implemented yet!');
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
                    xBounce(4:6) = obj.vel + dtBounce * obj.ballFlightModel3D(obj.vel);
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
                velNext = xBounce(4:6) + dt * obj.ballFlightModel3D(xBounce(4:6));
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
                    xContact(4:6) = obj.vel + dtContact * obj.ballFlightModel3D(obj.vel);
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
                velNext = xContact(4:6) + dt * obj.ballFlightModel3D(xContact(4:6));
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