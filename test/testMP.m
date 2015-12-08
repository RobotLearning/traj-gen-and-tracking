%% Testing mp for 2d-case

clc; clear; close all

%{
b2 = 2;
posRobotInit = [0;0];
posBallInit = [0;b2];
v1 = 1;
v2 = -1;
velBallInit = [v1;v2];
robotInit = [posRobotInit];
ballInit = [posBallInit;velBallInit];
solve_method = 'BVP';

% calculate ball path
dt = 0.1;
t = dt:dt:2;
N = length(t);
ballTime = t;
ballPath = posBallInit*ones(1,N) + velBallInit*t;
ballVel = velBallInit*ones(1,N);
ballPred = [ballPath;ballVel];

[t,y,u,J] = mp0(robotInit,ballTime,ballPred);

p1 = -v1;
speed = sqrt(v1^2 + v2^2);
T = b2 / speed

figure(1);
plot(y(1,:),y(2,:),'b');
hold on;
grid on;
plot(ballPath(1,:),ballPath(2,:),'r');
xlabel('x');
ylabel('y');
legend('robot','ball');
hold off;

%% Testing minimum principle for simple 3d- ball interception
 
% simple dynamics scenario to catch an incoming ball
posRobotInit = [0;0;0];
posBallInit = [0;0;3];
velRobotInit = [0;0;0];
velBallInit = [1;1;-1];
robotInit = [posRobotInit;velRobotInit];
ballInit = [posBallInit;velBallInit];
solve_method = 'BVP';

% calculate ball path
g = -9.8;
dt = 0.02;
t = dt:dt:1;
N = length(t);
ballTime = t;
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[0;0;g]*(t.^2);
ballVel = velBallInit*ones(1,N) + [0;0;g]*t;
ballPred = [ballPath; ballVel];
desVel = -ballVel;

[t,y,u,J] = mp(robotInit,ballTime,ballPred,desVel,solve_method);


figure(1);
plot3(y(1,:),y(2,:),y(3,:),'b');
hold on;
grid on;
plot3(ballPath(1,:),ballPath(2,:),ballPath(3,:),'r');
xlabel('x');
ylabel('y');
zlabel('z');
legend('robot','ball');
hold off;

figure(2);
plot(t,y(1:3,:)','-'); 
legend('x','y','z');
figure(3);
plot(t, u, '-');
legend('u1','u2','u3');
%axis([0 time(1,end) -1.5 3]);

% Solve symbolically to compare

% syms p1 p2 p3 p4 p5 p6 T v1 v2 v3 b1 b2 b3 
% eqns = [-1/6*p1*T^3 - 1/2*p4*T^2 == v1*T + b1, ... % desired ball pos at T
%         -1/6*p2*T^3 - 1/2*p5*T^2 == v2*T + b2, ...
%         -1/6*p3*T^3 - 1/2*p6*T^2 == v3*T + 1/2*g*T^2 + b3, ...
%         -1/2*p1*T^2 - p4*T == -v1, ... % desired racket velocities at T
%         -1/2*p2*T^2 - p5*T == -v2, ...
%         -1/2*p3*T^2 - p6*T == -v3 - g*T, ...
%          p1*v1 + p2*v2 + p3*(v3+g*T) - p6*g + ...
%          +1/2*(p4^2 + p5^2 + p6^2) + ...
%          (p1*v1 + p2*v2 + p3*(v3 + g*T)) == 0]; 
% S = solve(eqns,[p1 p2 p3 p4 p5 p6 T]);\
% T = subs(S.T,[v1,v2,v3,b1,b2,b3],[velBallInit',posBallInit']);
% T = eval(real(vpa(T)));
% fprintf('T = %f.\n', max(T));
%}

%% Testing MP on RR

clc; clear; close all;

% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 4;
% dimension of the output y
SIM.dimy = 4;
% dimension of the control input
SIM.dimu = 2;
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
l1 = 0.50; %length of first link, m
l2 = 0.40; %length of second link, m
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m2*l2^2; %kg m^2
% motor parameters
J_a1 = 0.100; % actuator inertia of link 1
J_g1 = 0.050; % gear inertia of link1
J_m1 = J_a1 + J_g1;
J_a2 = 0.080; % actuator inertia of link 2
J_g2 = 0.040; % gear inertia of link2
J_m2 = J_a2 + J_g2;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r_1 = 20;
r_2 = 20;
 
% pass it as a parameter structure
PAR.const.g = g;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;

% TODO: we should observe joint angles only
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON.link1.u.min = -Inf;
CON.link1.u.max = Inf;
CON.link2.u.min = -Inf;
CON.link2.u.max = Inf;
CON.link1.udot.min = -Inf;
CON.link1.udot.max = Inf;
CON.link2.udot.min = -Inf;
CON.link2.udot.max = Inf;

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,0,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
rr = RR(PAR,CON,COST,SIM);

% simple dynamics scenario to catch an incoming ball
robotInit = [0;0;0;0];
% ball is in 2d space
posBallInit = [1;0.5];
velBallInit = [-1;-1];
ballInit = [posBallInit;velBallInit];
%solve_method = 'BVP';

% calculate ball path
dt = 0.02;
t = dt:dt:1;
N = length(t);
ballTime = t;
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[g;0]*(t.^2);
ballVel = velBallInit*ones(1,N) + [g;0]*t;
ballPred = [ballPath; ballVel];
desVel = -ballVel;

% supply kinematics and jacobian
kinFnc = @(q) RRKinematics(q,PAR);

jacExact = @(q) [-l1*sin(q(1)) - l2*sin(q(1)+q(2)), -l2*sin(q(1)+q(2));
                 l1*cos(q(1)) + l2*cos(q(1)+q(2)), l2*cos(q(1)+q(2))];
derJacTimesQdot = @(q,qdot) [-l1*cos(q(1))*qdot(1) - l2*cos(q(1)+q(2))*(qdot(1)+qdot(2)), ...
                             -l2*cos(q(1)+q(2))*(qdot(1)+qdot(2));
                             -l1*sin(q(1))*qdot(1) - l2*sin(q(1)+q(2))*(qdot(1)+qdot(2)), ...
                             -l2*sin(q(1)+q(2))*(qdot(1)+qdot(2))];                             
%xdotExact = jacExact(q)*qdot;             
             
% test if jacobians match
% z = [0,0,1;
%      0,0,1];
% % get position of endeffector
% [j1,endEffPos] = kinFnc(q);
% endEffPos = [endEffPos;0];
% j1 = [j1;0];
% jointPos = [zeros(1,3);j1'];            
% 
% jac = jacobian(endEffPos',jointPos,z);
% xdot = jac(1:2,1:2)*qdot;


[t,y,u,J] = mpq(robotInit,ballTime,ballPred,desVel,jacExact,kinFnc,derJacTimesQdot);

%rr.animateArm(qact(1:2,:),ref);