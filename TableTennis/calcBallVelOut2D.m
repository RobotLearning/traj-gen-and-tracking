% Find the desired outgoing velocity of the ball
% At a certain ball position for 2D table tennis
% to land the ball at time2reach seconds to ballDes location (y-z)
%
% Reverting the nonlinear flight model and solving BVP to calculate 
% the initial outgoing velocity vout

function velOut = calcBallVelOut2D(ballDes,ballPos,time2reach,par)

    fast = par.fast;
    gravity = par.g;
    Cdrag = par.Cdrag;

    if fast
        % desired pos is in the centre of opponents court
        ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
        ballOutVel(2) = (ballDes(2) - ballPos(2) - ...
                        0.5*gravity*time2reach^2)/time2reach;
        % hack for now
        velOut = [1.1;1.2] .* ballOutVel(:);
        return;
    end

    % boundary value condition
    bc = @(x0,xf) [x0(1) - ballPos(1);
                   x0(2) - ballPos(2);
                   xf(1) - ballDes(1);
                   xf(2) - ballDes(2)];
    meshpoints = 50;
    solinit = bvpinit(linspace(0,time2reach,meshpoints),...
                      @(t)linFlightTraj(t,ballPos,ballDes,time2reach,par));
    sol = bvp4c(@(t,x)flightModel(t,x,par),bc,solinit);
    ballOut = deval(sol,0);
    velOut = ballOut(3:4);

end

% Flight model for integration
function xdot = flightModel(t,x,par)
    
    Cdrag = par.Cdrag;
    gravity = par.g;

    xdot(1) = x(3);
    xdot(2) = x(4);
    xdot(3) = -Cdrag * x(3) * sqrt(x(3)^2 + x(4)^2);
    xdot(4) = gravity - Cdrag * x(4) * sqrt(x(3)^2 + x(4)^2);

end

% initialize using a linear model (no drag)
function ball = linFlightTraj(t,ballPos,ballDes,time2reach,par)

    gravity = par.g;
    % desired pos is in the centre of opponents court
    ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
    ballOutVel(2) = (ballDes(2) - ballPos(2) - ...
                    0.5*gravity*time2reach^2)/time2reach;

    ball(1:2) = ballPos(1:2) + ballOutVel(:)*t;
    ball(3:4) = ballOutVel(:);

    ball(2) = ball(2) + 0.5*gravity*t^2;
    ball(4) = ball(4) + gravity*t;

end