%% Ball flight model and symplectic integration functions

% incorporate bounce also
function xNext = symplecticFlightModel(x,dt,params)

C = params.C;
g = params.g;
zTable = params.zTable;
yNet = params.yNet;
tableLength = params.table_length;
tableWidth = params.table_width;
% coeff of restitution-friction vector
K(1) = params.CFTX;
K(2) = params.CFTY;
K(3) = params.CRT;

xNext = zeros(6,1);
xNext(4:6) = x(4:6) + dt * ballFlightModel(x(4:6),C,g);
xNext(1:3) = x(1:3) + dt * xNext(4:6);

% condition for bouncing
if xNext(3) < zTable && abs(xNext(2) - yNet) < tableLength/2 && abs(xNext(1)) < tableWidth/2
    tol = 1e-4;
    dt1 = 0;
    dt2 = dt;
    xBounce = x;
    % doing bisection to find the bounce time
    while abs(xBounce(3) - zTable) > tol
        dtBounce = (dt1 + dt2) / 2;
        xBounce(4:6) = x(4:6) + dtBounce * ballFlightModel(x(4:6),C,g);
        xBounce(1:3) = x(1:3) + dtBounce * xBounce(4:6);
        if xBounce(3) > zTable
            % increase the time
            dt1 = dtBounce;
        else
            dt2 = dtBounce;
        end
    end
    % rebound
    xBounce(4:6) = reboundModel(xBounce(4:6),K);
    % integrate for the rest
    dt = dt - dtBounce;
    xNext(4:6) = xBounce(4:6) + dt * ballFlightModel(xBounce(4:6),C,g);
    xNext(1:3) = xBounce(1:3) + dt * xNext(4:6);
end

end

function xddot = ballFlightModel(xdot,C,g)

v = sqrt(xdot(1)^2 + xdot(2)^2 + xdot(3)^2);
xddot(1) = -C * v * xdot(1);
xddot(2) = -C * v * xdot(2);
xddot(3) = g - C * v * xdot(3);

xddot = xddot(:);
end

% K is the coefficient values in x-y-z directions
function xdot = reboundModel(xdot,K)

if K(3) > 0
    K(3) = -K(3); % make sure value is below zero
end

M = diag(K);
xdot = M * xdot;
    
end
