%% Testing minimum principle for simple ball interception
 
% simple dynamics scenario to catch an incoming ball
posRobotInit = [0;0;0];
posBallInit = [0;0;3];
velRobotInit = [0;0;0];
velBallInit = randn(3,1);
robotInit = [posRobotInit;velRobotInit];
ballInit = [posBallInit;velBallInit];
Tinit = 2.0;
solve_method = 'BVP';

[t,y,u,J] = mp(robotInit,ballInit,Tinit,solve_method);

% calculate ball path
g = -9.8;
N = length(t);
ball = posBallInit*ones(1,N) + velBallInit*t + 0.5*[0;0;g]*(t.^2);

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
plot(t,y(1:3,:)','-'); 
legend('x','y','z');
figure(3);
plot(t, u, '-');
legend('u1','u2','u3');
%axis([0 time(1,end) -1.5 3]);

%% Solve symbolically to compare

syms p1 p2 p3 p4 p5 p6 T v1 v2 v3 b1 b2 b3 
eqns = [-1/12*p1*T^3 - 1/4*p4*T^2 == v1*T + b1, ...
        -1/12*p2*T^3 - 1/4*p5*T^2 == v2*T + b2, ...
        -1/12*p3*T^3 - 1/4*p6*T^2 == v3*T + 1/2*g*T^2 + b3, ...
        -1/4*p1*T^2 - 1/2*p4*T == -v1, ...
        -1/4*p2*T^2 - 1/2*p5*T == -v2, ...
        -1/4*p3*T^2 - 1/2*p6*T == -v3 - g*T, ...
         p1*v1 + p2*v2 + p3*(v3+g*T) == 0]; 
S = solve(eqns,[p1 p2 p3 p4 p5 p6 T]);
T = subs(S.T,[v1,v2,v3,b1,b2,b3],[velBallInit',posBallInit']);
T = eval(real(vpa(T)));
fprintf('T = %f.\n', max(T));