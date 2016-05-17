%% Test table tennis functions here

clc; clear; close all;

loadTennisTableValues();

%% Testing the code for calculating outgoing velocities

% ballDes(1) = 0.0;
% ballDes(2) = dist_to_table - 3*table_y/2;
% ballDes(3) = table_z + ball_radius;

ballDes(1) = 0.0;
ballDes(2) = dist_to_table - table_length;
ballDes(3) = table_z + ball_radius + 0.2;

ballPos = [0;0;0];
time2reach = 0.5;
fast = false;

% calculate using BVP
par.fast = fast;
par.g = gravity;
par.Cdrag = Cdrag;
par.zTable = table_z;
par.radius = ball_radius;
par.bounce = true;
par.CFTX = CFTX; 
par.CFTY = CFTY; 
par.CRT = CRT;

ballOutVel = calcBallVelOut3D(ballDes,ballPos,time2reach,par);

% now check with symplectic Euler
dt = 0.02;
timeInt = 0;

ballFlightFnc = @(x) [x(4:6);ballFlightModel(x(4:6),Cdrag,gravity)];
x = [ballPos;ballOutVel(:)];

tSpan = [0,time2reach];
[T,X] = ode45(@(t,x)ballFlightFnc(x(:)),tSpan,x');
des_ode45 = X(end,1:3)

while timeInt < time2reach
    
    % Runge Kutta
    k1 = dt * ballFlightFnc(x);
    x_k1 = x + k1/2;
    k2 = dt * ballFlightFnc(x_k1);
    x_k2 = x + k2/2;
    k3 = dt * ballFlightFnc(x_k2);
    x_k3 = x + k3;
    k4 = dt * ballFlightFnc(x_k3);
    x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    
    % Symplectic Euler
    %{
    xdot = x(4:6);
    v = sqrt(xdot(1)^2 + xdot(2)^2 + xdot(3)^2);
    xddot(1) = -Cdrag * v * xdot(1);
    xddot(2) = -Cdrag * v * xdot(2);
    xddot(3) = gravity - Cdrag * v * xdot(3);
    xddot = xddot(:);
    
    x(4:6) = x(4:6) + dt * xddot(1:3);
    x(1:3) = x(1:3) + dt * x(4:6);
    %}
    
    timeInt = timeInt + dt;
end

% should be close to ballDes
des_symplecticEuler = x(1:3)