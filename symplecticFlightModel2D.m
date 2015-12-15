%% Ball flight model and symplectic integration functions
% Note: can also be run backwards in time by specifying -dt
% TODO: bounce model cannot be run backwards as of now!

% incorporate bounce also
function xNext = symplecticFlightModel2D(x,dt,params)

C = params.C;
g = params.g;
zTable = params.zTable;
yNet = params.yNet;
tableLength = params.table_length;
% coeff of restitution-friction vector
K(1) = params.CFTY;
K(2) = params.CRT;

xNext = zeros(4,1);
xNext(3:4) = x(3:4) + dt * ballFlightModel(x(3:4),C,g);
xNext(1:2) = x(1:2) + dt * xNext(3:4);

% condition for bouncing
if xNext(2) < zTable && abs(xNext(1) - yNet) < tableLength/2
    tol = 1e-4;
    dt1 = 0;
    dt2 = dt;
    xBounce = x;
    dtBounce = 0.0;
    % doing bisection to find the bounce time
    while abs(xBounce(2) - zTable) > tol
        dtBounce = (dt1 + dt2) / 2;
        xBounce(3:4) = x(3:4) + dtBounce * ballFlightModel(x(3:4),C,g);
        xBounce(1:2) = x(1:2) + dtBounce * xBounce(3:4);
        if xBounce(2) > zTable
            % increase the time
            dt1 = dtBounce;
        else
            dt2 = dtBounce;
        end
    end
    % rebound
    xBounce(3:4) = reboundModel(xBounce(3:4),K);
    % integrate for the rest
    dt = dt - dtBounce;
    xNext(3:4) = xBounce(3:4) + dt * ballFlightModel(xBounce(3:4),C,g);
    xNext(1:2) = xBounce(1:2) + dt * xNext(3:4);
end

end

function xddot = ballFlightModel(xdot,C,g)

v = sqrt(xdot(1)^2 + xdot(2)^2);
xddot(1) = -C * v * xdot(1);
xddot(2) = g - C * v * xdot(2);

xddot = xddot(:);
end

% K is the coefficient values in y-z directions
function xdot = reboundModel(xdot,K)

if K(2) > 0
    K(2) = -K(2); % make sure value is below zero
end

M = diag(K);
xdot = M * xdot;
    
end
