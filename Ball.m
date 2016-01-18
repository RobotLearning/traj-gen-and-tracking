%% Ball class for simulating all sorts of ball interactions

% 2d as well as 3d models

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
        % K is the coefficient values in y-z directions
        K
        % table structure containing table tennis-related values
        TABLE
        % reset event in case of landing, net contact, hit floor, etc.
        reset
    end
    
    % methods that all balls share
    methods (Access = public)   
        
        function obj = Ball(PAR)
            
            obj.C = PAR.C;
            obj.g = PAR.g;
            
            % table related values
            obj.TABLE.z = PAR.zTable;
            obj.TABLE.NET.y = PAR.yNet;
            obj.TABLE.length = PAR.table_length;
            % coeff of restitution-friction vector
            K(1) = PAR.CFTY;
            K(2) = PAR.CRT;
            
            if K(2) > 0
                K(2) = -K(2); % make sure value is below zero
            end
            
            obj.TABLE.K = K;
        end
        
        % 2D flight model
        function xddot = ballFlightModel2D(obj)

            v = sqrt(obj.vel(1)^2 + obj.vel(2)^2);
            xddot(1) = -obj.C * v * obj.vel(1);
            xddot(2) = obj.g - obj.C * v * obj.vel(2);

            xddot = xddot(:);
        end
        
        % Symplectic Integration for 2D case
        function symplecticFlightModel2d(obj,dt)            

            vell = obj.vel + dt * obj.ballFlightModel2D();
            poss = obj.pos + dt * vell;
             
            if obj.checkNet(dt)
                % ball will touch the net
                % TODO
            end
            if obj.checkRebound(dt)
                % get rebound time                
                % TODO
            end
            if obj.checkContactWithRacket(dt)
            else
                obj.vel = vell;
                obj.pos = poss;
            end
            
        end
        
        % K is the coefficient values in y-z directions
        function reboundModel(obj)

            M = diag(obj.TABLE.K);
            obj.vel = M * obj.vel;
        end
        
        % check if rebound event happens within dt seconds
        function checkRebound(obj,dt)
            % condition for bouncing

            % adding tolerance to make sure nothing weird happens when the ball
            % is below the table
            d = 0.1; 

            % if ball is above the table and on the table plane
            % consider rebound only
            if obj.pos(2) < obj.TABLE.z && obj.pos(2) > obj.TABLE.z - d && ...
               abs(obj.pos(1) - obj.TABLE.NET.y) < obj.TABLE.length/2
                tol = 1e-4;
                dt1 = 0;
                dt2 = dt;
                xBounce = obj.pos;
                dtBounce = 0.0;
                % doing bisection to find the bounce time
                while abs(xBounce(2) - obj.TABLE.z) > tol
                    dtBounce = (dt1 + dt2) / 2;
                    xBounce(3:4) = x(3:4) + dtBounce * ballFlightModel(x(3:4),C,g);
                    xBounce(1:2) = x(1:2) + dtBounce * xBounce(3:4);
                    if xBounce(2) > zTable
                        % increase the time
                        dt1 = dtBounce;
                    else
                        dt2 = dtBounce;
                    end
                end
                % rebound
                xBounce(3:4) = reboundModel(xBounce(3:4),K);
                % integrate for the rest
                dt = dt - dtBounce;
                xNext(3:4) = xBounce(3:4) + dt * ballFlightModel(xBounce(3:4),C,g);
                xNext(1:2) = xBounce(1:2) + dt * xNext(3:4);
            end
        end
        
    end
    
end