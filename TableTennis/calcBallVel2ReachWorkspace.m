% Find the desired outgoing velocity of the ball
% At a certain ball position for 3D table tennis
% to reach desired workspace location
%
% Reverting the nonlinear flight model and solving BVP to calculate 
% the initial outgoing velocity vout
%
% TODO: not working! Rewrite without using bvp's!

function velOut = calcBallVel2ReachWorkspace(ballDes,ballPos,time2reach,par)

    fast = par.fast;
    gravity = par.g;
    Cdrag = par.Cdrag;

    if fast
        % desired pos is in the centre of opponents court
        ballOutVel(1) = (ballDes(1) - ballPos(1))/time2reach;
        ballOutVel(2) = (ballDes(2) - ballPos(2))/time2reach;
        ballOutVel(3) = (ballDes(3) - ballPos(3) - ...
                        0.5*gravity*time2reach^2)/time2reach;
        % hack for now
        velOut = [1.1;1.1;1.2] .* ballOutVel(:);
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
                      @(t)linBounceTraj(t,ballPos,ballDes,time2reach,par));
    sol = bvp4c(@(t,x)flightModel(t,x,par),bc,solinit);
    ballOut = deval(sol,0);
    velOut = ballOut(4:6);

end

% Flight model for integration
function xdot = flightModel(t,x,par)
    
    Cdrag = par.Cdrag;
    gravity = par.g;
    zTable = par.zTable;
    ballRadius = par.radius;

    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = -Cdrag * x(4) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
    xdot(5) = -Cdrag * x(5) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
    xdot(6) = gravity - Cdrag * x(6) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);

    if x(3) < zTable + ballRadius && xdot(3) < 0.0 && ...
            x(2) < par.dist_to_table
        M = diag(par.coeff_bounce);
        velPreBounce = xdot(1:3);
        xdot(1:3) = M * velPreBounce(:);
    end

end

% initialize using bounce + linear model (no drag)
function ball = linBounceTraj(t,ballPos,ballDes,time2reach,par)

    gravity = par.g;
    zTable = par.zTable;
    ballRadius = par.radius;
    M = par.coeff_bounce;

    table_z = zTable + ballRadius;
    a = 0.5*gravity*(1 - M(3));
    b = time2reach * gravity * (0.5*M(3) - 1);
    c = 0.5*gravity*time2reach^2 - M(3)*(table_z-ballPos(3)) + (table_z - ballDes(3));
    d = M(3)*time2reach*(table_z-ballPos(3));
    Ts = roots([a b c d]);
    % get the positive ones
    tland_est = Ts(imag(Ts) == 0 & (Ts > 0.2));
    tland_est = max(tland_est);
    
    ballOutVel(3) = (table_z - ballPos(3) - 0.5*gravity*tland_est^2)/tland_est;
    ballOutVel(2) = (ballDes(2) - ballPos(2)) / (tland_est + M(2)*(time2reach-tland_est));
    ballOutVel(1) = (ballDes(1) - ballPos(1)) / (tland_est + M(1)*(time2reach-tland_est));
    
    if t < tland_est
        ball(1:3) = ballPos(1:3) + ballOutVel(:)*t;
        ball(4:6) = ballOutVel(:);
        ball(3) = ball(3) + 0.5*gravity*t^2;
        ball(6) = ball(6) + gravity*t;
    else
        % find landing x-y locations
        xt = ballPos(1) + tland_est*ballOutVel(1);
        yt = ballPos(2) + tland_est*ballOutVel(2);
        ball(1) = xt + (t - tland_est)*M(1)*ballOutVel(1);
        ball(2) = yt + (t - tland_est)*M(2)*ballOutVel(2);
        ball(3) = table_z + (t - tland_est)*M(3)*(ballOutVel(3) + gravity*tland_est);
        ball(4:6) = M .* ballOutVel(:);
        ball(3) = ball(3) + 0.5*gravity*(t - tland_est)^2;
        ball(6) = ball(6) + gravity*(t - tland_est);
    end

end