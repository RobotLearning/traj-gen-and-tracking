%% TODO: Testing mp for 2d-case

clc; clear; close all

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

%% Solve symbolically to compare

syms p1 p2 p3 p4 p5 p6 T v1 v2 v3 b1 b2 b3 
eqns = [-1/6*p1*T^3 - 1/2*p4*T^2 == v1*T + b1, ... % desired ball pos at T
        -1/6*p2*T^3 - 1/2*p5*T^2 == v2*T + b2, ...
        -1/6*p3*T^3 - 1/2*p6*T^2 == v3*T + 1/2*g*T^2 + b3, ...
        -1/2*p1*T^2 - p4*T == -v1, ... % desired racket velocities at T
        -1/2*p2*T^2 - p5*T == -v2, ...
        -1/2*p3*T^2 - p6*T == -v3 - g*T, ...
         p1*v1 + p2*v2 + p3*(v3+g*T) - p6*g + ...
         +1/2*(p4^2 + p5^2 + p6^2) + ...
         (p1*v1 + p2*v2 + p3*(v3 + g*T)) == 0]; 
S = solve(eqns,[p1 p2 p3 p4 p5 p6 T]);
T = subs(S.T,[v1,v2,v3,b1,b2,b3],[velBallInit',posBallInit']);
T = eval(real(vpa(T)));
fprintf('T = %f.\n', max(T));