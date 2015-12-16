%% Testing mp for 2d-case

clc; clear; close all

%{
b1 = 0;
b2 = 2;
x1 = 0; x2 = 0;
posRobotInit = [x1;x2];
posBallInit = [b1;b2];
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

% speed = sqrt(v1^2 + v2^2);
% T = b2 / speed;

x0 = [1;1;0.5];
x = fsolve(@(x) bv1(x(1),x(2),x(3),x1,x2,b1,b2,v1,v2),x0);
T = x(end)

syms p1 p2 T;
eqns = [-p1*T + x1 == b1 + v1*T, ...
        -p2*T + x2 == b2 + v2*T, ...
        1/2*(p1^2 + p2^2) + p1*v1 + p2*v2 == 0];
S = solve(eqns,[p1,p2,T]);
T = subs(S.T);
T = eval(real(vpa(T)));
fprintf('T = %f.\n', max(T));

figure(1);
plot(y(1,:),y(2,:),'b');
hold on;
grid on;
plot(ballPath(1,:),ballPath(2,:),'r');
xlabel('x');
ylabel('y');
legend('robot','ball');
hold off;
%}

%% Testing minimum principle for simple 3d- ball interception
 
% simple dynamics scenario to catch an incoming ball
%{
posRobotInit = [0;0;0];
b1 = 0; b2 = 0; b3 = 3;
posBallInit = [b1;b2;b3];
v1 = 1; v2 = 1; v3 = -1;
velRobotInit = [0;0;0];
velBallInit = [v1;v2;v3];
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
syms p1 p2 p3 p4 p5 p6 T 
eqns = [1/6*p1*T^3 + 1/2*p4*T^2 == v1*T + b1, ... % desired ball pos at T
        1/6*p2*T^3 + 1/2*p5*T^2 == v2*T + b2, ...
        1/6*p3*T^3 + 1/2*p6*T^2 == v3*T + 1/2*g*T^2 + b3, ...
        1/2*p1*T^2 + p4*T == -v1, ... % desired racket velocities at T
        1/2*p2*T^2 + p5*T == -v2, ...
        1/2*p3*T^2 + p6*T == -v3 - g*T, ...
         p1*v1 + p2*v2 + p3*(v3+g*T) + (p6+p3*T)*g + ...
         +1/2*(p4^2 + p5^2 + p6^2) == 0]; 
S = solve(eqns,[p1 p2 p3 p4 p5 p6 T]);
T = subs(S.T);
T = eval(real(vpa(T)));
fprintf('Symbolic solver: T = %f.\n', max(T));

x0 = [zeros(6,1);0.5];
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
    'MaxFunEvals',1000,'MaxIter',1000);
lb = [-Inf(6,1);0];
x = lsqnonlin(@(x) bv1(x,[b1;b2;b3;v1;v2;v3]), x0, lb, [], options);
T = x(end)

%}

%% Testing MP on RR

%{
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
g = 1; %9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
l1 = 1.0; %length of first link, m
l2 = 1.0; %length of second link, m
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
q10 = 0.0; q20 = 0.0;
q0 = [q10;q20];
robotInit = [q0;0;0];
% ball is in 2d space
b1 = 1; b2 = 1;
v1 = -0.2; v2 = -0.2;
posBallInit = [b1; b2];
velBallInit = [v1; v2];
ballInit = [posBallInit;velBallInit];
%solve_method = 'BVP';

% calculate ball path
dt = 0.02;
t = dt:dt:1;
N = length(t);
ballTime = t;
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[g;0]*(t.^2);
ballVel = velBallInit*ones(1,N);
ballPred = [ballPath; ballVel];
desVel = -ballVel;
% 
% derJacTimesQdot = @(q,qdot) [-l1*cos(q(1))*qdot(1) - l2*cos(q(1)+q(2))*(qdot(1)+qdot(2)), ...
%                              -l2*cos(q(1)+q(2))*(qdot(1)+qdot(2));
%                              -l1*sin(q(1))*qdot(1) - l2*sin(q(1)+q(2))*(qdot(1)+qdot(2)), ...
%                              -l2*sin(q(1)+q(2))*(qdot(1)+qdot(2))];                             
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


% solve numerically with lsqnonlin
x0 = [pi/4*ones(4,1);0.5];
lb = [-Inf(4,1);0];

PAR.ball.x0 = [b1;b2];
PAR.ball.v0 = [v1;v2];
PAR.ball.g = g; % acceleration
PAR.robot.Q0 = q0;
PAR.robot.Qd0 = 0;
PAR.robot.class = rr;

tic;
if 1 %isempty(which('lsqnonlin'))
    % try newton raphson
    options = optimset('TolX',1e-12); % set TolX
    fun = @(x) bvMP(x,PAR);
    [x, resnorm, resval, exitflag, output, jacob] = newtonraphson(fun, x0, options);
    fprintf('\nExitflag: %d, %s\n',exitflag, output.message) % display output message
else
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
        'MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-20);
    [x,resnorm,resval] = lsqnonlin(@(x) bvMP(x,PAR), x0, lb, [], options);
end
toc
T = x(end)

t = dt:dt:T;
N = length(t);
mu0 = x(1:2);
nu0 = -x(3:4) - mu0*T;
q = (1/6)*mu0*t.^3 + (1/2)*nu0*t.^2 + repmat(q0,1,N);
%[~,xc] = kinFnc(q);
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[g;0]*(t.^2);
%plot(x(2,:),-x(1,:),'b',ballPath(2,:),-ballPath(1,:),'r');

%[t,y,u,J] = mpq(robotInit,ballTime,ballPred,desVel,jacExact,kinFnc,derJacTimesQdot);
%t(end)

rr.animateArm(q(1:2,:),ballPath);

% Solve symbolically to compare
% syms p1 p2 p3 p4 T 
% eqns = [l1*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) + ... % kinematics constr.
%         l2*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) == b1 + v1*T, ...
%         l1*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) + ... % kinematics constr.
%         l2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) == b2 + v2*T, ...
%         (-l1*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) - ... % Jacobian constr.
%         l2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)) * ((-1/2)*p1*T^2 - p3*T) + ...
%         (l2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)) * ((-1/2)*p3*T^2 - p4*T) == -v1, ...
%          (l1*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) + ... % Jacobian constr.
%         l2*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)) * ((-1/2)*p1*T^2 - p3*T) + ...
%         (l2*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)) * ((-1/2)*p3*T^2 - p4*T) == -v2, ...
%          -1/2*(p3^2 + p4^2) == (1/(l1*l2*(cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1)* ... %Hamiltonian
%          sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) - ...
%          sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) * ...
%                -sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)))) ...
%                * l2*v1*p1*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) + p1*l2*v2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) - p2*l1*v1*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) ...
%                -p2*l2*v1*cos(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2) - p2*l1*v2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1) ...
%                -p2*l2*v2*sin(-(1/6)*p1*T^3 - (1/2)*p3*T^2 + q1 ...
%                -(1/6)*p2*T^3 - (1/2)*p4*T^2 + q2)];
% S = solve(eqns,[p1 p2 p3 p4 T]);
% T = subs(S.T);
% T = eval(real(vpa(T)));
% fprintf('T = %f.\n', max(T));

%}

%% Testing MP for RRR for hitting an uncertain ball

%%{
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
g = 1; %9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
m3 = 0.5; 
l1 = 1.0; %length of first link, m
l2 = 1.0; %length of second link, m
l3 = 1.0;
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

% simple dynamics scenario to catch an incoming ball
q10 = 0.0; q20 = 0.0; q30 = 0.0;
q0 = [q10;q20;q30];
robotInit = [q0;0;0;0];
% ball is in 2d space
b1 = 1; b2 = 1;
v1 = -0.2; v2 = -0.2;
posBallInit = [b1; b2];
velBallInit = [v1; v2];
ballInit = [posBallInit;velBallInit];
%solve_method = 'BVP';

% give ball path as a bunch of points or as a function of x0,v0 and T
dt = 0.02;
Tmax = 5.0;
t = dt:dt:Tmax;
N = length(t);
ballTime = t;
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[g;0]*(t.^2);
ballVel = velBallInit*ones(1,N);
ballPred = [ballPath; ballVel];
desVel = -ballVel;
% ball path as a function
ballFnc = @(b0,v0,T) [b0 + v0*T + (1/2)*[g;0]*T^2; v0 + [g;0]*T];

% supply kinematics and jacobian
% kinFnc = @(q) RRRKinematics(q,PAR);
% 
% jacExact = @(q) [-l1*sin(q(1)) - l2*sin(q(1)+q(2)) - l3*sin(q(1)+q(2)+q(3)), ...
%                  -l2*sin(q(1)+q(2)) - l3*sin(q(1)+q(2)+q(3)), ...
%                  -l3*sin(q(1)+q(2)+q(3)); 
%                  l1*cos(q(1)) + l2*cos(q(1)+q(2)) + l3*cos(q(1)+q(2)+q(3)), ...
%                  l2*cos(q(1)+q(2)) + l3*cos(q(1)+q(2)+q(3)), ...
%                  l3*cos(q(1)+q(2)+q(3))];
% derJacTimesQdot = @(q,qdot) [-l1*cos(q(1))*qdot(1) - l2*cos(q(1)+q(2))*(qdot(1)+qdot(2)), ...
%                              -l2*cos(q(1)+q(2))*(qdot(1)+qdot(2));
%                              -l1*sin(q(1))*qdot(1) - l2*sin(q(1)+q(2))*(qdot(1)+qdot(2)), ...
%                              -l2*sin(q(1)+q(2))*(qdot(1)+qdot(2))];                                 


% solve numerically with lsqnonlin
x0 = [pi/4*ones(6,1);0.5];
lb = [-Inf(6,1);0];

PAR.ball.x0 = [b1;b2];
PAR.ball.v0 = [v1;v2];
PAR.ball.g = g; % acceleration
PAR.ball.model = ballFnc;
PAR.ball.time = ballTime;
PAR.ball.path = ballPred;
PAR.robot.Q0 = q0;
PAR.robot.Qd0 = 0;
PAR.robot.class = rrr;

tic;
%options = optimset('TolX',1e-12); % set TolX
%fun = @(x) bvMP(x,PAR);
%[x, resnorm, resval, exitflag, output, jacob] = newtonraphson(fun, x0, options);
%fprintf('\nExitflag: %d, %s\n',exitflag, output.message) % display output message
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
     'MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-20);
[x,resnorm,resval] = lsqnonlin(@(x) bvMP(x,PAR), x0, lb, [], options);
toc
T = x(end)

t = dt:dt:T;
N = length(t);
mu0 = x(1:3);
nu0 = -x(4:6) - mu0*T;
q = (1/6)*mu0*t.^3 + (1/2)*nu0*t.^2 + repmat(q0,1,N);
%[~,xc] = kinFnc(q);
N = length(t);
ballPath = posBallInit*ones(1,N) + velBallInit*t + 0.5*[g;0]*(t.^2);

rrr.animateArm(q(1:3,:),ballPath);
%}