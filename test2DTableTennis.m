%% Test 2D table tennis with RRR arm
% There won't be any X-coordinates in the 2d version
% First considering hitting the ball only
clc; clear; close all;

% load table parameters
loadTennisTableValues;
ball(1:2,1) = ball_cannon(2:3);
ball(3:4,1) = [4.000 3.2] + 0.05 * randn(1,2);
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
l1 = 0.4; %length of first link, m
l2 = 0.4; %length of second link, m
l3 = 0.4;
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
eps = 1e-1;
C = [eye(2),zeros(2)];

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
% coeff of restitution-friction vector
params.CFTY = CFTY;
params.CRT = CRT;

funState = @(x,u,dt) symplecticFlightModel2D(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(2);
filter = EKF(dim,funState,mats);

%% Prepare the animation

% simple dynamics scenario to catch an incoming ball
q10 = 0.3; q20 = -0.3; q30 = -0.4;
q0 = [q10;q20;q30];
% get joints, endeffector pos and orientation
[x1,x2,x3,mats] = rrr.kinematics(q0);

% rotate everything by 90 degrees
R = [0 1; -1 0];
x1 = R * x1;
x2 = R * x2;
x3 = R * x3;
shift = [0;0];
link1_y = [shift(1) x1(1,1)];
link1_z = [shift(2) x1(2,1)];
link2_y = [x1(1,1),x2(1,1)];
link2_z = [x1(2,1),x2(2,1)];
link3_y = [x2(1,1),x3(1,1)];
link3_z = [x2(2,1),x3(2,1)];
hf = figure('color','white');
axis equal; axis auto;
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
h7 = line(line_racket_y,line_racket_z,'color',[.4 .4 .4],'LineWidth',4);
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
predict = false;
% maximum horizon to predict
time2PassTable = 1.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;

% initialize the filters state with sensible values
ballNoisyVel = ball(3:4,1);% + sqrt(eps)*rand(2,1);
filter.initState([ball(1:2,1);ballNoisyVel],eps);

% initialize indices and time
i = 1; k = 1; t = 0.0;
dt = 0.01;
ballNoisyPos = ball(1:2,1);

while numTrials < 50
    
    % Evolve ball flight, table and racket contact interactions
    if ball(2,i) <= floor_level
        disp('Ball hit the floor. Resetting...');
        t = 0.0; i = 1; k = 1;
        predict = false;
        numTrials = numTrials+1;
        time2PassTable = 1.0;
        ball(1:2,1) = ball_cannon(2:3);
        ball(3:4,1) = [4.000 3.2] + 0.05 * randn(1,2);
        ballNoisyVel = ball(3:4,1);% + sqrt(eps)*rand(2,1);
        filter.initState([ball(1:2,1);ballNoisyVel],eps);
    end
    
    ball(:,i+1) = funState(ball(:,i),0,dt);
    
    %TODO: check contact with racket
    tol = ball_radius;
    l = min(k,size(path,2));
    vecFromRacketToBall = ball(1:2,i) - robot(5:6,l);
    racketPlane = racket_dir(:,l);
    projPlane = racketPlane*racketPlane'/(racketPlane'*racketPlane);
    projOrth = eye(2) - projPlane;
    distToRacketPlane = norm(projOrth * vecFromRacketToBall);
    distOnRacketPlane = norm(projPlane * vecFromRacketToBall);
    if distToRacketPlane < tol && distOnRacketPlane < racket_radius
        disp('A hit! Well done!');
        numHits = numHits+1;
        % change ball velocity based on contact model
    end
    
    % Get noisy ball positions till ball passes the net
    if time2PassTable >= maxTime2Hit
        % add noise
        ballNoisyPos = ball(1:2,i+1); % + sqrt(eps) * randn(2,1);
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
        % predict ball trajectory once
        if ~predict            
            predictHorizon = maxPredictHorizon;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(4,predictLen);
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
            ballIncoming = ballPred; %ballPred(:,idxAfterTable);
            ballTime = (1:predictLen) * dt; %idxAfterTable * dt;
            minTimeToHit = ballTime(1);
            
            % COMPUTE OPTIMAL TRAJECTORY HERE
            % first computing in cartesian space assuming a cartesian robot
            T0 = 0.5;
            B0 = ballIncoming(:,1);
            
            % solve numerically with Newton-Raphson / nlsq-optimization
            x0 = [pi/4*ones(6,1);T0];
            lb = [-Inf(6,1);0];

            PAR.ball.x0 = B0(1:2);
            PAR.ball.v0 = B0(3:4);
            PAR.ball.g = g; % acceleration
            %PAR.ball.model = @(b,v,T) funState([b;v],0,T);
            PAR.ball.time = ballTime;
            PAR.ball.path = ballPred;
            PAR.robot.Q0 = q0;
            PAR.robot.Qd0 = 0;
            PAR.robot.class = rrr;
            PAR.rotate = R;
            
            tic;
            %{
            options = optimset('TolX',1e-12); % set TolX
            fun = @(x) bvMP(x,PAR);
            [x, resnorm, resval, exitflag, output, jacob] = newtonraphson(fun, x0, options);
            fprintf('\nExitflag: %d, %s\n',exitflag, output.message) % display output message
            %}
            %%{
            options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
                 'MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-20);
            [x,resnorm,resval] = lsqnonlin(@(x) bvMP(x,PAR), x0, lb, [], options);
            %fprintf('Residuals:\n');
            %fprintf('%f\n', resval); 
            %}
            toc
            T = x(end)

            t = dt:dt:T;
            N = length(t);
            mu0 = x(1:3);
            nu0 = -x(4:6) - mu0*T;
            q = (1/6)*mu0*t.^3 + (1/2)*nu0*t.^2 + repmat(q0,1,N);
            qd = (1/2)*mu0*t.^2 + nu0*t;
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

        % Initiate the robot
        robot(:,k+1) = path(:,min(k+1,size(path,2)));
        k = k + 1;
        
    end 
    
    % ANIMATE BOTH THE ROBOT AND THE BALL
    set(h0,'XData',ball(1,i));
    set(h0,'YData',ball(2,i));    
    
    
    link1_y = [shift(1),robot(1,k)];
    link1_z = [shift(2),robot(2,k)];
    link2_y = [robot(1,k),robot(3,k)];
    link2_z = [robot(2,k),robot(4,k)];
    link3_y = [robot(3,k),robot(5,k)];
    link3_z = [robot(4,k),robot(6,k)];
    line_racket_y = [robot(7,k),robot(9,k)];
    line_racket_z = [robot(8,k),robot(10,k)];
     
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
    i = i + 1;
    
    drawnow;
    pause(0.002);
        
end

hold off;