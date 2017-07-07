%% Ball spin model and symplectic integration functions

function x_next = discr_const_spin_model(x,dt,params)

w0 = params.init_spin;
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

x_next = zeros(6,1);

switch alg
    case 'Euler'
        x_next(4:6) = x(4:6) + dt * ball_const_spin_model(x(4:6),Cdrag,Clift,w0,g);
        x_next(1:3) = x(1:3) + dt * x_next(4:6);
    case 'RK4'
        ballFlightFnc = @(x) [x(4:6);ball_const_spin_model(x(4:6),Cdrag,Clift,w0,g)];
        k1 = dt * ballFlightFnc(x);
        x_k1 = x + k1/2;
        k2 = dt * ballFlightFnc(x_k1);
        x_k2 = x + k2/2;
        k3 = dt * ballFlightFnc(x_k2);
        x_k3 = x + k3;
        k4 = dt * ballFlightFnc(x_k3);
        x_next = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    otherwise
        error('Not implemented!');
end

% condition for bouncing
if x_next(3) < zTable + ballRadius && ...
        abs(x_next(2) - yNet) < tableLength/2 && abs(x_next(1)) < tableWidth/2 ...
        && x_next(6) < 0.0
    tol = 1e-4;
    dt1 = 0;
    dt2 = dt;
    x_bounce = x;
    dt_bounce = 0.0;
    iter = 0;
    % doing bisection to find the bounce time
    while iter < 5 %abs(xBounce(3) - zTable) > tol
        dt_bounce = (dt1 + dt2) / 2;
        x_bounce(4:6) = x(4:6) + dt_bounce * ball_const_spin_model(x(4:6),Cdrag,Clift,w0,g);
        x_bounce(1:3) = x(1:3) + dt_bounce * x_bounce(4:6);
        if x_bounce(3) > zTable
            % increase the time
            dt1 = dt_bounce;
        else
            dt2 = dt_bounce;
        end
        iter = iter + 1;
    end
    % rebound
    x_bounce(4:6) = rebound_const_spin_model(x_bounce(4:6),w0,e_t,mu,ballRadius);
    % integrate for the rest
    dt = dt - dt_bounce;
    x_next(4:6) = x_bounce(4:6) + dt * ball_const_spin_model(x_bounce(4:6),Cdrag,Clift,w0,g);
    x_next(1:3) = x_bounce(1:3) + dt * x_next(4:6);
end

end

% K is the coefficient values in x-y-z directions
function xdot = rebound_const_spin_model(xdot,w0,e_t,mu,r)

    vbT = [xdot(1) - r*w0(2);
           xdot(2) + r*w0(1);
           0];

    %e_t = 0.90;
%     nu_s = 1 - (2/5)*mu*(1+e_t)*abs(xdot(3))/norm(vbT);
%     alpha = 2/5; % roll
%     if nu_s > 0 % slide 
%         alpha = mu * (1 + e_t) * abs(xdot(3))/norm(vbT);
%     end
    % REPLACING ABOVE SINCE BALL ALWAYS SLIDES IT SEEMS
    alpha = mu * (1 + e_t) * abs(xdot(3))/norm(vbT)

    Av = diag([1.0-alpha,1.0-alpha,-e_t]); 
    Bv = [0, alpha*r, 0; 
          -alpha*r, 0, 0; 
          zeros(1,3)];
    Aw = [0, -3*alpha/(2*r), 0;
          3*alpha/(2*r), 0, 0;
          zeros(1,3)];
    Bw = diag([1-3*alpha/2, 1-3*alpha/2, 1]);

    %M = [Av, Bv; Aw, Bw];
    M = [Av,Bv];
    xdot = M * [xdot;w0(:)];
    
end