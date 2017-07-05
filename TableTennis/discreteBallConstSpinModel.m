%% Ball spin model and symplectic integration functions

% spinNext is modified if there is a bounce
function [xNext,spinNext] = discreteBallConstSpinModel(x,dt,params)

spin = params.spin;
Clift = params.Clift;
Cdrag = params.C;
g = params.g;
zTable = params.zTable;
yNet = params.yNet;
ballRadius = params.radius;
tableLength = params.table_length;
tableWidth = params.table_width;
e_t = params.CRT; 
mu = params.mu; % dynamic coeff of friction
alg = params.ALG;

% now we also include spin
xNext = zeros(6,1);

switch alg
    case 'Euler'
        acc = ballSpinModel([x(4:6);spin(:)],Cdrag,Clift,g);
        xNext(4:6) = x(4:6) + dt * acc(1:3);
        xNext(1:3) = x(1:3) + dt * xNext(4:6);
    case 'RK4'
        ballFlightFnc = @(x) [x(4:6);ballSpinModel([x(4:6);spin(:)],...
                                      Cdrag,Clift,g)];
        k1 = dt * ballFlightFnc(x);
        x_k1 = x + k1(1:6)/2;
        k2 = dt * ballFlightFnc(x_k1);
        x_k2 = x + k2(1:6)/2;
        k3 = dt * ballFlightFnc(x_k2);
        x_k3 = x + k3(1:6);
        k4 = dt * ballFlightFnc(x_k3);
        xNext = x + (k1(1:6) + 2*k2(1:6) + 2*k3(1:6) + k4(1:6))/6;
    otherwise
        error('Not implemented!');
end

% condition for bouncing
if xNext(3) < zTable + ballRadius && ...
        abs(xNext(2) - yNet) < tableLength/2 && abs(xNext(1)) < tableWidth/2
    tol = 1e-4;
    dt1 = 0;
    dt2 = dt;
    xBounce = x;
    dtBounce = 0.0;
    iter = 0;
    % doing bisection to find the bounce time
    while iter < 5 %abs(xBounce(3) - zTable) > tol
        dtBounce = (dt1 + dt2) / 2;
        acc = ballSpinModel([xBounce(4:6);spin(:)],Cdrag,Clift,g);
        xBounce(4:6) = x(4:6) + dtBounce * acc(1:3);
        xBounce(1:3) = x(1:3) + dtBounce * xBounce(4:6);
        if xBounce(3) > zTable
            % increase the time
            dt1 = dtBounce;
        else
            dt2 = dtBounce;
        end
        iter = iter + 1;
    end
    % rebound
    xdot = reboundSpinModel([xBounce(4:6);spin(:)],e_t,mu,ballRadius);
    xBounce(4:6) = xdot(1:3);
    spinNext = xdot(4:6);
    % integrate for the rest
    dt = dt - dtBounce;
    acc = ballSpinModel([xBounce(4:6);spin(:)],Cdrag,Clift,g);
    xNext(4:6) = xBounce(4:6) + dt * acc(1:3);
    xNext(1:3) = xBounce(1:3) + dt * xNext(4:6);
end

end