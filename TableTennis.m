%% Table Tennis class for simulating a match between two robots!

% TODO: this should be part of a TASK superclass?
classdef TableTennis < handle
    
    properties
        
        % ball class
        ball
        % robot 1
        robot1
        % robot 2 (opponent)
        robot2
        
    end
    
    methods
        
        %% CONSTRUCTOR
        function obj = TableTennis(wam,wam2)
            
            % initialize two robots
            obj.robot1 = wam;
            obj.robot2 = wam2;
            % initialize a ball
            obj.ball = Ball();             
            % initialize animation
            obj.initAnimation(q0,qd0);
            
        end
        
        %% Animation functions here
        function initAnimation(obj,q0,qd0)           

            % Prepare the animation
            loadTennisTableValues();
            wam = obj.robot1;
            b = obj.ball;

            % get joints, endeffector pos and orientation
            [x,xd,o] = wam.calcRacketState([q0;qd0]);
            [joint,ee,racket] = wam.drawPosture(q0);
            q = q0; qd = qd0;

            figure;
            %uisetcolor is useful here
            orange = [0.9100 0.4100 0.1700];
            gray = [0.5020    0.5020    0.5020];
            ballSurfX = b.pos(1) + b.MESH.X;
            ballSurfY = b.pos(2) + b.MESH.Y;
            ballSurfZ = b.pos(3) + b.MESH.Z;
            % transform into base coord.
            h1 = surf(ballSurfX,ballSurfY,ballSurfZ);
            set(h1,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
            hold on;
            h10 = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
            endeff = [joint(end,:); ee];
            h11 = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
            h12 = fill3(racket(1,:), racket(2,:), racket(3,:), 'r');

            %h2 = scatter3(ballPred(1,1),ballPred(2,1),ballPred(3,1),20,'b','filled');
            %h3 = scatter3(X0(1,1),X0(2,1),X0(3,1),20,'k','filled');
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
            net_width = 0.01;
            
        end
        
        % update the animation for robot and the ball
        function updateAnimation(obj)
            
            % ANIMATE BOTH THE ROBOT AND THE BALL
            b = obj.ball;
    
            % ball
            set(h1,'XData',b.pos(1) + b.MESH.X);
            set(h1,'YData',b.pos(2) + b.MESH.Y);
            set(h1,'ZData',b.pos(3) + b.MESH.Z);

            % robot joints
            set(h10,'XData',joint(:,1));
            set(h10,'YData',joint(:,2));
            set(h10,'ZData',joint(:,3));
            % robot endeffector
            set(h11,'XData',endeff(:,1));
            set(h11,'YData',endeff(:,2));
            set(h11,'ZData',endeff(:,3));
            % robot racket
            set(h12,'XData',racket(1,:));
            set(h12,'YData',racket(2,:));
            set(h12,'ZData',racket(3,:));

            drawnow;
            pause(0.001);
            
        end
        
        %% MAIN LOOP
        
        % a robot serves for 2 points
        % and then they switch
        % till one gets 11 points or more with 2 difference
        function match(obj)
            
            b = obj.ball;
            r1 = obj.robot1;
            r2 = obj.robot2;
            % points for the robots
            points1 = 0;
            points2 = 0;
            
            while max(points1,points) < 11 || abs(points1-points2) < 2
                
                ballObs = obj.getBallObs();
                obj.play(r1,r2,ballObs);
                obj.judge(points1,points2);
                
                
            end
            
        end
        
        function play(obj,robot,robot2,ballObs)
           
            % if ball is NaN serve the ball
            
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end