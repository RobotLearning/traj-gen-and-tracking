%% Ball spin model and symplectic integration functions

function xNext = discreteBallSpinModel(x,dt,params)

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
xNext = zeros(12,1);

switch alg
    case 'Euler'
        xNext(7:12) = x(7:12) + dt * ballSpinModel(x(7:12),Cdrag,Clift,g);
        xNext(1:6) = x(1:6) + dt * xNext(7:12);
    case 'RK4'
        ballFlightFnc = @(x) [x(7:12);ballSpinModel(x(7:12),Cdrag,Clift,g)];
        k1 = dt * ballFlightFnc(x);
        x_k1 = x + k1/2;
        k2 = dt * ballFlightFnc(x_k1);
        x_k2 = x + k2/2;
        k3 = dt * ballFlightFnc(x_k2);
        x_k3 = x + k3;
        k4 = dt * ballFlightFnc(x_k3);
        xNext = x + (k1 + 2*k2 + 2*k3 + k4)/6;
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
        xBounce(7:12) = x(7:12) + dtBounce * ballSpinModel(x(7:12),Cdrag,Clift,g);
        xBounce(1:6) = x(1:6) + dtBounce * xBounce(7:12);
        if xBounce(3) > zTable
            % increase the time
            dt1 = dtBounce;
        else
            dt2 = dtBounce;
        end
        iter = iter + 1;
    end
    % rebound
    xBounce(7:12) = reboundSpinModel(xBounce(7:12),e_t,mu,ballRadius);
    % integrate for the rest
    dt = dt - dtBounce;
    xNext(7:12) = xBounce(7:12) + dt * ballSpinModel(xBounce(7:12),Cdrag,Clift,g);
    xNext(1:6) = xBounce(1:6) + dt * xNext(7:12);
end

end

% K is the coefficient values in x-y-z directions
function xdot = reboundSpinModel(xdot,e_t,mu,r)

    vbT = [xdot(1) - r*xdot(5);
           xdot(2) + r*xdot(4);
           0];

    nu_s = 1 - (2/5)*mu*(1+e_t)*abs(xdot(3))/norm(vbT);
    alpha = 2/5; % roll
    if nu_s > 0 % slide 
        alpha = mu * (1 + e_t) * abs(xdot(3))/norm(vbT);
    end

    Av = diag([1.2,1-alpha,-e_t]); % 1-alpha for x also in the paper!!!
    Bv = [0, alpha*r, 0;
          -alpha*r, 0, 0;
          zeros(1,3)];
    Aw = [0, -3*alpha/(2*r), 0;
          3*alpha/(2*r), 0, 0;
          zeros(1,3)];
    Bw = diag([1-3*alpha/2, 1-3*alpha/2, 1]);

    M = [Av, Bv; Aw, Bw];
    xdot = M * xdot;
    %M = [Av,Bv];
    %xdot(1:3) = M*xdot;

    
end
