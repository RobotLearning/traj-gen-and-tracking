% Find the desired outgoing velocity of the ball
% At a certain predicted ball pos for 2D table tennis
%
% Reverting the nonlinear flight model and solving BVP to calculate 
% the initial outgoing velocity vout

function vout = calcBallVelOut2D(ballDes,ballPos,time2reach)

loadTennisTableValues();

% approximate with a linear model
ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
ballOutVel(2) = (ballDes(2) - ballPos(2) - ...
                    0.5*gravity*time2reach^2)/time2reach;

% initialize using a linear model (no drag)
linFlightTraj = @(t) [ballPos + ballOutVel(:)*t;
                      ballOutVel(:)] + ...
                     [0;0.5*gravity*t^2;0;gravity*t];
flightModel = @(t,x) [x(3);
                      x(4);                                  
                     -Cdrag * x(3) * sqrt(x(3)^2 + x(4)^2);
                     gravity - Cdrag * x(4) * sqrt(x(3)^2 + x(4)^2)];
% boundary value condition
bc = @(x0,xf) [x0(1) - ballPos(1);
               x0(2) - ballPos(2);
               xf(1) - ballDes(1);
               xf(2) - ballDes(2)];
meshpoints = 50;
solinit = bvpinit(linspace(0,time2reach,meshpoints),...
                   linFlightTraj);
sol = bvp4c(flightModel,bc,solinit);
ballOut = deval(sol,0);
vout = ballOut(3:4);