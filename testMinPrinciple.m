%% Testing minimum principle for table tennis
 
function [T,J] = testMinPrinciple()
 
% simple dynamics scenario to catch an incoming ball
posRobotInit = [0;0;0];
posBallInit = [0;0;3];
velRobotInit = [0;0;0];
velBallInit = randn(3,1);
momentumInit = zeros(6,1);

% try to get also an analytic soln.
% extend to t0 variable as well!
tic
Tinit = 2.0;
solinit = bvpinit(linspace(0,1),[posRobotInit;velRobotInit;
                                 momentumInit;Tinit]);

bcfull = @(y0,yf) bc(y0,yf,posBallInit,velBallInit,...
                     posRobotInit,velRobotInit);                             
sol = bvp4c(@ode, bcfull, solinit);
y = sol.y;
time = y(end)*sol.x;
% Calculate u from lambda
lambda = y(7:12,:);
B = [zeros(3);eye(3)];
R = eye(3);
u = -R \ (B' * lambda);
toc

% calculate ball path
g = -9.8;
N = length(time);
ball = posBallInit*ones(1,N) + velBallInit*time + 0.5*[0;0;g]*(time.^2);

% Calculate cost
T = time(end);
dt = T/N;
J = 0.5 * dt * trace(u'*R*u);

figure(1);
plot3(y(1,:),y(2,:),y(3,:),'b');
hold on;
grid on;
plot3(ball(1,:),ball(2,:),ball(3,:),'r');
xlabel('x');
ylabel('y');
zlabel('z');
legend('robot','ball');
hold off;

figure(2);
plot(time,y(1:3,:)','-'); 
legend('x','y','z');
figure(3);
plot(time, u, '-');
legend('u1','u2','u3');
%axis([0 time(1,end) -1.5 3]);

end

% -------------------------------------------------------------------------
% ODE's of augmented states
function ydot = ode(t,y)
    
    T = y(end);
    x = y(1:6);
    p = y(7:12);
    A = [zeros(3),eye(3); zeros(3,6)];
    B = [zeros(3);eye(3)];
    R = 1;
    u = -R \ (B' * p);
    x_dot = T*(A*x + B*u);
    p_dot = -A' * p;
    ydot = T*[x_dot; p_dot; 0];
end

% -------------------------------------------------------------------------
% boundary conditions: 
%
function res = bc(y0,yf,posBallInit,velBallInit,posRobotInit,velRobotInit)
    % ball at time T
    T = yf(end);
    g = -9.8;
    bT = posBallInit + velBallInit*T + 0.5*T^2*[0;0;g];
    dbdT = velBallInit + T*[0;0;g];
    % 6 initial conditions
    res = [ y0(1) - posRobotInit(1);
            y0(2) - posRobotInit(2);
            y0(3) - posRobotInit(3);
            y0(4) - velRobotInit(1);
            y0(5) - velRobotInit(2);
            y0(6) - velRobotInit(3);
            yf(1:3) - bT; 
            yf(4:6) + dbdT;
            %yf(10:12) - 0;
            yf(7:9)' * dbdT];
end

% % 1. Estimate initial pos and vel with kalman smoother up to net/bounce
% % 2. Predict the ball trajectory after bounce
% % 3. Fix the center of the table as the goal position
% %    Solve BVP problems (bvp4c) by sampling from the predicted traj.
% % 4. Using the contact model estimate racket vels. and feasible
% % orientations and interpolate through whole traj.
% %
% % 5. Using an optimal control toolbox
% % State the initial position of the robot
% % State the input penalties
% % State the transversality conditions
% %
% % Transversality conditions:
% % - positions of the ball should touch the racket (tubular region around
% %   the path)
% % - velocities of the racket are fixed
% % - orientation of the racket could be fixed (although there is some
% % freedom)
% % - only using kinematics and jacobian
% %
% % 6. Compare with existing virtual hitting plane method
% % - construct VHP
% % - estimate pos and vel. as well as desired vel.
% % - use inverse kinematics and polynomials to construct traj.
% %
% % 7. Evolve both with the nominal model (inverse dynamics) and compare
% %   probability of returning the ball given actual ball velocities 
% %   and/or control errors
% 
% %% Load table values
% 
% % load table parameters
% loadTennisTableValues;
% 
% %% Initialize Barrett WAM
% 
% initializeWAM;
% 
% %% initialize EKF
% dim = 6;
% eps = 1e-3;
% C = [eye(3),zeros(3)];
% 
% params.C = Cdrag;
% params.g = gravity;
% params.zTable = table_z;
% params.yNet = dist_to_table - table_y;
% params.table_length = table_length;
% params.table_width = table_width;
% % coeff of restitution-friction vector [TO BE LEARNED!]
% params.CFTX = CFTX;
% params.CFTY = CFTY;
% params.CRT = CRT;
% 
% funState = @(x,u,dt) symplecticFlightModel(x,dt,params);
% % very small but nonzero value for numerical stability
% mats.O = eps * eye(dim);
% mats.C = C;
% mats.M = eps * eye(3);
% filter = EKF(dim,funState,mats);
% 
% %% Main control and estimation loop
% 
% % flags for the main loop
% finished = false;
% predict = false;
% % maximum horizon to predict
% maxPredictHorizon = 0.8;
% % land the ball on the centre of opponents court
% desBall(1) = 0.0;
% desBall(2) = dist_to_table - 3*table_y/2;
% desBall(3) = table_z + ball_radius;
% time2reach = 0.5;
% % initialize ball on the ball cannon with a sensible vel
% ball(1:3,1) = ball_cannon;
% ball(4:6,1) = [-0.9 4.000 3.2] + 0.05 * randn(1,3);
% % initialize the filters state with sensible values
% filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);
% 
% % initialize indices and time
% i = 1;
% dt = 0.01;
% ballNoisyPos = ball(1:3,1);
% % time to simulate in total
% timeSimTotal = 1.2;
% t = 0.0;
% 
% while t < timeSimTotal && ~finished
%     
%     ball(:,i+1) = funState(ball(:,i),0,dt);
%     
%     %% Get noisy ball positions till ball passes the net
%     if ballNoisyPos(2) <= dist_to_table - table_y - ball_radius;
%         % add noise
%         ballNoisyPos = ball(1:3,i+1) + sqrt(eps) * randn(3,1);
%         filter.linearize(dt,0);
%         filter.update(ballNoisyPos,0);
%         filter.predict(dt,0);
%         % robot just waits here
%     else
%         if ~predict
%             %% predict ball trajectory
%             % ideally stop predicting only when time left is
%             % less than a certain value, e.g. 0.5 sec
%             predictHorizon = min(maxPredictHorizon,timeSimTotal - t);
%             predictLen = floor(predictHorizon / dt);
%             ballPred = zeros(6,predictLen);
%             for j = 1:predictLen
%                 filter.linearize(dt,0);
%                 filter.predict(dt,0);
%                 ballPred(:,j) = filter.x;
%             end
%             predict = true;
%         end
%         
%         %% Calculate the intersection of ball path with robot workspace
%         
%         % for now only considering the ball positions after table
%         tol = 5e-2;
%         idxAfterTable = find(ballPred(2,:) > dist_to_table + tol);
%         ballIncoming = ballPred(:,idxAfterTable);
%         timeIncoming = idxAfterTable * dt;
%         minTimeToHit = timeIncoming(1);
%         ballTrjPred = [timeIncoming; ballIncoming];
%         
%         %% Calculate ball outgoing velocities attached to each ball pos
%         
%         for j = 1:size(ballTrjPred,2)
%             % desired pos is in the centre of opponents court
%             velBallOut(1) = (desBall(1) - ballIncoming(1,j))/time2reach;
%             velBallOut(2) = (desBall(2) - ballIncoming(2,j))/time2reach;
%             velBallOut(3) = (desBall(3) - ballIncoming(3,j) - ...
%                             0.5*gravity*time2reach^2)/time2reach;
%             % initialize using a linear model (no drag)
%             linFlightTraj = @(t) [ballIncoming(1:3,j) + velBallOut(:)*t;
%                                   velBallOut(:)] + ...
%                                   [0;0;0.5*gravity*t^2;0;0;gravity*t];
%             flightModel = @(t,x) [x(4);
%                                 x(5);
%                                 x(6);
%                           -Cdrag * x(4) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
%                           -Cdrag * x(5) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
%                   gravity - Cdrag * x(6) * sqrt(x(4)^2 + x(5)^2 + x(6)^2)];
%             % boundary value condition
%             bc = @(x0,xf) [x0(1) - ballIncoming(1,j);
%                            x0(2) - ballIncoming(2,j);
%                            x0(3) - ballIncoming(3,j);
%                            xf(1) - desBall(1);
%                            xf(2) - desBall(2);
%                            xf(3) - desBall(3)];
%             meshpoints = 50;
%             solinit = bvpinit(linspace(0,time2reach,meshpoints),...
%                               linFlightTraj);
%             sol = bvp4c(flightModel,bc,solinit);
%             ballOut = deval(sol,0);
%             velOut(:,j) = ballOut(4:6);
%         
%             %% Use the inverse contact model to compute racket vels and normal
%             % at every point
%             
%             % FOR NOW USING ONLY THE MIRROR LAW
%             normal(:,j) = (velOut(:,j) - ballIncoming(4:6,j)) ...
%                           ./ norm(velOut(:,j) - ballIncoming(4:6,j),2);
%             velOutAlongNormal = velOut(:,j)' * normal(:,j);
%             velInAlongNormal = ballIncoming(4:6,j)' * normal(:,j);
%             racketVelAlongNormal(j) = velOutAlongNormal + ... 
%                 CRR * velInAlongNormal / (1 + CRR);
%         end
%         
%         
%         %% OPTIMAL CONTROL HERE
%         
%         
%     end
%             
%         
%         
%     t = t + dt;
%     i = i + 1;
% end
% 
% %% Plot the balls path
% 
% figure(4);
% scatter3(ball(1,:),ball(2,:),ball(3,:),'r');
% hold on;
% %scatter3(r1,r2,r3,'b');
% scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
% %scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
% title('Ball-robot interaction');
% grid on;
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% legend('ball', 'ballPred');
% %legend('ball','robot');
% fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
% fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
% hold off;
