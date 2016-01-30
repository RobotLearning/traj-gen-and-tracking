%% Test 2D table tennis with RRR arm
% There won't be any X-coordinates in the 2d version
clc; clear; close all;

% load table parameters
loadTennisTableValues;
ball(1:2,1) = ball_cannon(2:3);
ball(3:4,1) = 1.1*[4.000 3.2] + 0.00 * randn(1,2);
ballPred(1:4,1) = ball(1:4,1); % to initialize drawing of predicted balls

%% Initialize the RRR robot
% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 6;
% dimension of the output y
SIM.dimy = 3;
% dimension of the control input
SIM.dimu = 3;
% time step h 
SIM.h = 0.01;
% measurement noise covariance
SIM.eps_m = 3e-10;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

% constants
g = 9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
m3 = 0.5; 
l1 = 0.5; %length of first link, m
l2 = 0.5; %length of second link, m
l3 = 0.5;
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
l_c3 = 0.5;
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m2*l2^2; %kg m^2
I3 = (1/12)*m3*l3^2;
% motor parameters
J_a1 = 0.100; % actuator inertia of link 1
J_g1 = 0.050; % gear inertia of link1
J_m1 = J_a1 + J_g1;
J_a2 = 0.080; % actuator inertia of link 2
J_g2 = 0.040; % gear inertia of link2
J_m2 = J_a2 + J_g2;
J_a3 = 0.12;
J_g3 = 0.06;
J_m3 = J_a3 + J_g3;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r_1 = 20;
r_2 = 20;
r_3 = 40;
 
% pass it as a parameter structure
PAR.const.g = g;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link3.mass = m3;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link3.length = l3;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link3.centre.dist = l_c3;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link3.inertia = I3;
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link3.motor.inertia = J_m3;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;
PAR.link3.motor.gear_ratio = r_3;

PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON.link1.u.min = -Inf;
CON.link1.u.max = Inf;
CON.link2.u.min = -Inf;
CON.link2.u.max = Inf;
CON.link3.u.min = -Inf;
CON.link3.u.max = Inf;
CON.link1.udot.min = -Inf;
CON.link1.udot.max = Inf;
CON.link2.udot.min = -Inf;
CON.link2.udot.max = Inf;
CON.link3.udot.min = -Inf;
CON.link3.udot.max = Inf;

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,1,0,0,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
rrr = RRR(PAR,CON,COST,SIM);

%% initialize EKF
dim = 4;
eps = 0.0001; %1e-1;
C = [eye(2),zeros(2)];

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
% coeff of restitution-friction vector
params.CFTY = CFTY;
params.CRT = CRT;

funState = @(x,u,dt) discreteBallFlightModel2D(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(2);
filter = EKF(dim,funState,mats);

%% Prepare the animation

% simple dynamics scenario to catch an incoming ball
q10 = 0.3; q20 = -0.3; q30 = -0.4;
q0 = [q10;q20;q30];
q(:,1) = q0;
qd(:,1) = zeros(3,1);
% get joints, endeffector pos and orientation
[x1,x2,x3,mats] = rrr.kinematics(q0);
x0 = x3(:,1);

% rotate everything by 90 degrees
R = [0 1; -1 0];
x1 = R * x1;
x2 = R * x2;
x3 = R * x3;
x0 = R * x0;
shift = [0;0];
link1_y = [shift(1) x1(1,1)];
link1_z = [shift(2) x1(2,1)];
link2_y = [x1(1,1),x2(1,1)];
link2_z = [x1(2,1),x2(2,1)];
link3_y = [x2(1,1),x3(1,1)];
link3_z = [x2(2,1),x3(2,1)];
hf = figure('color','white');
axis equal; %axis auto;
%uisetcolor is useful here
orange = [0.9100 0.4100 0.1700];
gray = [0.5020    0.5020    0.5020];
h0 = scatter(ball(1,1),ball(2,1),20,orange,'filled');
hold on;
h1 = line(link1_y, link1_z, 'color', [.4 .4 .4],'LineWidth',4);
h2 = line(link2_y, link2_z, 'color', [.4 .4 .4],'LineWidth',4);
h3 = line(link3_y, link3_z, 'color', [.4 .4 .4],'LineWidth',4);
h4 = scatter(shift(1),shift(2),50,'k','LineWidth',4);
h5 = scatter(x1(1,1),x1(2,1),20,'b','LineWidth',4);
h6 = scatter(x2(1,1),x2(2,1),20,'b','LineWidth',4);
%h7 = scatter(x3(1,1),x3(2,1),20,'b','LineWidth',4);
% add the endeffector as a racket-line attached
racket_orient = R * mats(1:2,1,3);
racket_dir = R * mats(1:2,2,3);
racket_radius = 0.08;
racket1 = x3 - racket_radius * racket_dir;
racket2 = x3 + racket_radius * racket_dir;
line_racket_y = [racket1(1,1),racket2(1,1)];
line_racket_z = [racket1(2,1),racket2(2,1)];
h7 = line(line_racket_y,line_racket_z,'color',[1 0 0],'LineWidth',4);
robot(:,1) = [x1;x2;x3;racket1;racket2];

%h3 = scatter3(X0(1,1),X0(2,1),X0(3,1),20,'k','filled');
title('Ball-robot interaction');
grid on;
axis equal;
xlabel('y');
ylabel('z');

% define table for 2d case
table_z = floor_level - table_height;
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

% tol_x = 0.1; tol_y = 0.1; tol_z = 0.3;
% xlim([dist_to_table - table_length - tol_y, tol_y]);
% ylim([table_z - tol_z, table_z + 2*tol_z]);
legend('ball','robot');
fill(T(:,1),T(:,2),[0 0.7 0.3]);
fill(net(:,1),net(:,2),[0 0 0]);

%% Main control and estimation loop

% flags for the main loop
contact = false;
numTrials = 0;
numHits = 0;
numLands = 0;
predict = false;
hit = false;
land = false;
net = false;
% maximum horizon to predict
time2PassTable = 1.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 1.0;
maxWaitAfterHit = 1.0;

% initialize the filters state with sensible values
ballNoisyVel = 1.1*[4.000; 3.2];
filter.initState([ball(1:2,1);ballNoisyVel],eps);

% initialize indices and time
ballIdx = 1; robotIdx = 1; t = 0.0;
dt = 0.01;
cnt = 0.0; % counter for resetting
ballNoisyPos = ball(1:2,1);

while numTrials < 50
    
    % Reset criterion
    if ball(2,ballIdx) <= floor_level || land || net %cnt > maxWaitAfterHit
        numTrials = numTrials + 1;
        fprintf('Success: %d/%d\n',numLands,numTrials);
        disp('Resetting...'); % ball hit the floor
        t = 0.0; ballIdx = 1; robotIdx = 1; cnt = 0.0;
        predict = false;
        hit = false;
        land = false;
        net = false;        
        time2PassTable = 1.0;
        ball(1:2,1) = ball_cannon(2:3);
        ball(3:4,1) = 1.1*[4.000; 3.2] + 0.0 * randn(2,1);
        ballNoisyVel = 1.1*[4.000; 3.2];
        filter.initState([ball(1:2,1);ballNoisyVel],eps);
    end

    % check contact with racket    
    tol = ball_radius;
    l = min(robotIdx,size(path,2));
    % get racket pos, vel and orientation (on plane)
    [xRacket,velRacket,mats] = rrr.getEndEffectorState(q(:,l),qd(:,l));
    xRacket = R * xRacket;
    velRacket = R * velRacket;
    racketPlane = R * squeeze(mats(1:2,2,3,:)); % orientation along racket
    vecFromRacketToBall = ball(1:2,ballIdx) - xRacket;
    %racketPlane = racket_dir(:,l);
    projPlane = racketPlane*racketPlane'/(racketPlane'*racketPlane);
    projOrth = eye(2) - projPlane;
    distOnRacketPlane = abs(racketPlane'*vecFromRacketToBall);
    distToRacketPlane = sqrt(abs(vecFromRacketToBall'*vecFromRacketToBall - distOnRacketPlane^2));
    if distToRacketPlane < tol && distOnRacketPlane < racket_radius && ~hit            
        %disp('A hit! Well done!');
        hit = true;
        fprintf('Hit at y = %f, z = %f\n',ball(1,ballIdx),ball(2,ballIdx));
        numHits = numHits + 1;
        % Change ball velocity based on contact model
        % get ball velocity
        velIn = ball(3:4,ballIdx);
        speedInAlongRacket = racketPlane'*velIn;
        velInAlongRacket = speedInAlongRacket'* racketPlane;
        velInAlongNormal = velIn - velInAlongRacket;
        speedRacketAlongRacket = racketPlane'*velRacket;
        velRacketAlongRacket = speedRacketAlongRacket'*racketPlane;
        velRacketAlongNormal = velRacket - velRacketAlongRacket;
        % this is kept the same in mirror law        
        velOutAlongNormal = velRacketAlongNormal + ...
            CRR * (velRacketAlongNormal - velInAlongNormal);
        velOut = velOutAlongNormal + velInAlongRacket;
        ball(3:4,ballIdx) = velOut;
    end    
    
    % run counter to terminate game after hit
    if hit
        cnt = cnt + dt;
        tol = 2*net_width_2d;
        % check contact with net
        if abs(ball(1,ballIdx) - (dist_to_table - table_y)) <= tol 
            if ball(2,ballIdx) < (table_z + net_height)
                disp('Hit the net! Resetting...');
                net = true;
            else
                disp('Passing over the net');
            end
        end
        % check for landing
        tol = 0.01;
        if ball(2,ballIdx) <= table_z + ball_radius + tol && ball(4,ballIdx) < 0             
            land = true;
            fprintf('Land at y = %f\n', ball(1,ballIdx));
            if abs(ball(1,ballIdx) - (dist_to_table - 3*table_y/2)) < table_y/2
                %disp('Ball landed! Amazing!'); 
                numLands = numLands + 1;                
                %else
                %disp('Ball is not inside the court. Lost a point!')
            end
        end
    end
    
    % evolve ball according to flight model
    ball(:,ballIdx+1) = funState(ball(:,ballIdx),0,dt);
    
    % Get noisy ball positions till ball passes the net
    if time2PassTable >= maxTime2Hit
        % add noise
        ballNoisyPos = ball(1:2,ballIdx+1); % + sqrt(eps) * randn(2,1);
        filter.linearize(dt,0);
        filter.update(ballNoisyPos,0);
        filter.predict(dt,0);
        xSave = filter.x;
        PSave = filter.P;
        % update the time it takes to pass table
        yBallEst = filter.x(1);
        tPredIncrement = dt;
        time2PassTable = 0.0;
        while yBallEst <= dist_to_table
             %filter.linearize(tPredIncrement,0);
             filter.predict(tPredIncrement,0);
             yBallEst = filter.x(1);
             time2PassTable = time2PassTable + tPredIncrement;
        end
        % revert back to saved state
        filter.initState(xSave,PSave);
    else
        % Predict ball trajectory once
        if ~predict            
            predictHorizon = maxPredictHorizon;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(4,predictLen);
            % FOR DEBUGGING
            filter.initState(ball(:,ballIdx),PSave);
            for j = 1:predictLen
                %filter.linearize(dt,0);
                filter.predict(dt,0);
                ballPred(:,j) = filter.x;
            end
            predict = true;

            % Calculate the intersection of ball path with robot workspace

            % for now only considering the ball positions after table
            tol = 5e-2;
            idxAfterTable = find(ballPred(1,:) > dist_to_table + tol);
            ballTime = (1:predictLen) * dt; %idxAfterTable * dt;
            ballPredWithTime = [ballPred;ballTime]; %ballPred(:,idxAfterTable);
            minTimeToHit = ballTime(1);
            
            %% COMPUTE TRAJECTORY HERE
            
            % consider different vhp trajectories and keep the optimal
            %{
            ballInitVar = PSave;
            [q,qd] = rrr.generateOptimalTTT(ballPred,ballTime,ballInitVar,desBall,time2reach,time2return,q0);
            %}
            
            %%{
            % Compute VHP Trajectory here            
            [q,qd,~] = rrr.generate2DTTTwithVHP(ballPred,ballTime,q0);            
%             M = 100;
%             ballInitVar = PSave;
%             bSamp = repmat(ballPred(:,1),1,M) + 0; %chol(ballInitVar)*randn(4,M);
%             tic;
%             probLand = calcProbOfLand(rrr,bSamp,q,qd)
%             toc
            %computeMPfor2DTT;            
            %}
            
            [x1,x2,x3,mats] = rrr.kinematics(q);
            x1 = R * x1;
            x2 = R * x2;
            x3 = R * x3;
            racket_dir = R * squeeze(mats(1:2,2,3,:));

            racket1 = x3 - racket_radius * racket_dir;
            racket2 = x3 + racket_radius * racket_dir;
            orient = [racket1;racket2];
            path = [x1;x2;x3;orient];
            
        end % end predict        

        % Robot joint coordinates in cartesian space (for drawing)
        robot(:,robotIdx+1) = path(:,min(robotIdx+1,size(path,2)));
        robotIdx = robotIdx + 1;
        
    end 
    
    %% ANIMATE BOTH THE ROBOT AND THE BALL
    set(h0,'XData',ball(1,ballIdx));
    set(h0,'YData',ball(2,ballIdx));    
    
    
    link1_y = [shift(1),robot(1,robotIdx)];
    link1_z = [shift(2),robot(2,robotIdx)];
    link2_y = [robot(1,robotIdx),robot(3,robotIdx)];
    link2_z = [robot(2,robotIdx),robot(4,robotIdx)];
    link3_y = [robot(3,robotIdx),robot(5,robotIdx)];
    link3_z = [robot(4,robotIdx),robot(6,robotIdx)];
    line_racket_y = [robot(7,robotIdx),robot(9,robotIdx)];
    line_racket_z = [robot(8,robotIdx),robot(10,robotIdx)];
     
    set(h1,'XData',link1_y)
    set(h1,'YData',link1_z)
    set(h2,'XData',link2_y)
    set(h2,'YData',link2_z)
    set(h3,'XData',link3_y)
    set(h3,'YData',link3_z)
    set(h5,'XData',link1_y(2));
    set(h5,'YData',link1_z(2));
    set(h6,'XData',link2_y(2));
    set(h6,'YData',link2_z(2));
    %set(h7,'XData',link3_y(2));
    %set(h7,'YData',link3_z(2));
    set(h7,'XData',line_racket_y);
    set(h7,'YData',line_racket_z);
    
    t = t + dt;
    ballIdx = ballIdx + 1;
    
    drawnow;
    pause(0.002);
        
end

hold off;