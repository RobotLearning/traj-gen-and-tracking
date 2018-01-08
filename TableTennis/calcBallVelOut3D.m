% Find the desired outgoing velocity of the ball
% At a certain ball position for 3D table tennis
% to land the ball at time2reach seconds to ballDes location (x-y-z)
%
% Reverting the nonlinear flight model and solving BVP to calculate 
% the initial outgoing velocity vout

function velOut = calcBallVelOut3D(ballDes,ballPos,time2reach,par)

    fast = par.fast;
    gravity = par.g;
    Cdrag = par.Cdrag;
    dt = 0.02;
    %t = dt:dt:time2reach;
    
    if fast
        % desired pos is in the centre of opponents court
        ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
        ballOutVel(2) = (ballDes(2) - ballPos(2))/time2reach;
        ballOutVel(3) = (ballDes(3) - ballPos(3) - ...
                        0.5*gravity*time2reach^2)/time2reach;
        % hack for now
        velOut = [1.1;1.1;1.2] .* ballOutVel(:);
        %ballOut = ballPos + velOut * t;
        %ballOut(3,:) = ballOut(3,:) + gravity * t/2;
        return;
    end

    % boundary value condition
    bc = @(x0,xf) [x0(1) - ballPos(1);
                   x0(2) - ballPos(2);
                   x0(3) - ballPos(3);
                   xf(1) - ballDes(1);
                   xf(2) - ballDes(2);
                   xf(3) - ballDes(3)];
    meshpoints = 50;
    solinit = bvpinit(linspace(0,time2reach,meshpoints),...
                      @(t)linFlightTraj(t,ballPos,ballDes,time2reach,par));
    sol = bvp4c(@(t,x)flightModel(t,x,par),bc,solinit);
    ballOut = deval(sol,0);
    velOut = ballOut(4:6);

end

% Flight model for integration
function xdot = flightModel(t,x,par)
    
    Cdrag = par.Cdrag;
    gravity = par.g;

    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = -Cdrag * x(4) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
    xdot(5) = -Cdrag * x(5) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
    xdot(6) = gravity - Cdrag * x(6) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);

end

% initialize using a linear model (no drag)
function ball = linFlightTraj(t,ballPos,ballDes,time2reach,par)

    gravity = par.g;
    % desired pos is in the centre of opponents court
    ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
    ballOutVel(2) = (ballDes(2) - ballPos(2))/time2reach;
    ballOutVel(3) = (ballDes(3) - ballPos(3) - ...
                    0.5*gravity*time2reach^2)/time2reach;

    ball(1:3) = ballPos(1:3) + ballOutVel(:)*t;
    ball(4:6) = ballOutVel(:);

    ball(3) = ball(3) + 0.5*gravity*t^2;
    ball(6) = ball(6) + gravity*t;

end