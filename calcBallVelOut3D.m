% Find the desired outgoing velocity of the ball
% At a certain ball position for 3D table tennis
%
% Reverting the nonlinear flight model and solving BVP to calculate 
% the initial outgoing velocity vout

function velOut = calcBallVelOut3D(ballDes,ballPos,time2reach,fast)

loadTennisTableValues();

% desired pos is in the centre of opponents court
ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
ballOutVel(2) = (ballDes(2) - ballPos(2))/time2reach;
ballOutVel(3) = (ballDes(3) - ballPos(3) - ...
                0.5*gravity*time2reach^2)/time2reach;
            
if fast
    velOut = ballOutVel;
    return;
end
            
% initialize using a linear model (no drag)
linFlightTraj = @(t) [ballPos(1:3) + ballOutVel(:)*t;
                      ballOutVel(:)] + ...
                      [0;0;0.5*gravity*t^2;0;0;gravity*t];
flightModel = @(t,x) [x(4);
                      x(5);
                      x(6);
                     -Cdrag * x(4) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
                     -Cdrag * x(5) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
                     gravity - Cdrag * x(6) * sqrt(x(4)^2 + x(5)^2 + x(6)^2)];
% boundary value condition
bc = @(x0,xf) [x0(1) - ballPos(1);
               x0(2) - ballPos(2);
               x0(3) - ballPos(3);
               xf(1) - ballDes(1);
               xf(2) - ballDes(2);
               xf(3) - ballDes(3)];
meshpoints = 50;
solinit = bvpinit(linspace(0,time2reach,meshpoints),linFlightTraj);
sol = bvp4c(flightModel,bc,solinit);
ballOut = deval(sol,0);
velOut = ballOut(4:6);