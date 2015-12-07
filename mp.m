%% Solving Minimum principle for 3d kinematics case
% Transversality conditions are applied
% TODO: test gradient descent also!

function [t,x,u,J] = mp(robotInit,ballTime,ballPred,desVel,solve_method)

    Tinit = 0.5;
    R = 1;
    momentumInit = zeros(6,1);
    tic;
    if strcmp(solve_method,'BVP')
        % Initial guess for the solution
        solinit = bvpinit(linspace(0,1),[robotInit;momentumInit;Tinit]);
        options = bvpset('Stats','on','RelTol',1e-1);

        bcfull = @(y0,yf) bc(y0,yf,ballTime,ballPred,robotInit,desVel);        
        sol = bvp4c(@ode, bcfull, solinit, options);
        y = sol.y;
        t = y(end)*sol.x;
        y = sol.y;
        x = y(1:6,:);
        % Calculate u from lambda
        p = y(7:12,:);
        B = [zeros(3);eye(3)];
        u = -R \ (B' * p);
        % Calculate cost
        T = t(end);
        N = length(t);
        dt = T/N;
        J = 0.5 * dt * trace(u'*R*u);
    else if strcmp(solve_method,'GD') % gradient descent    
        %% Steepest descent
        max_iteration = 100;
        N = 100;
        J = zeros(1,max_iteration);
        eps = 1e-3;
        step = 0.2;
        options = [];
        tu = linspace(0,1,N);
        u = zeros(1,N);
        for i = 1:max_iteration
            % 1) start with assumed control u and move forward
            initx = [posRobotInit;velRobotInit;Tinit];
            [t,x] = ode45(@(t,x) stateEq(t,x,u,tu), [0 1], ...
                           initx, options);
            x = x(1:6,:);
            % calculate initp
            initp = [momentumInit;0];
            [tp,p] = ode45(@(tp,p) costateEq(tp,p,u,tu,x,t), [0 1], initp, options);
            % Important: costate is stored in reverse order. The dimension of
            % costates may also different from dimension of states
            % Use interploate to make sure x and p is aligned along the time axis
            p = interp1(tp,p,t);
            % Calculate deltaH 
            dH = Hu(x1,p,t,u,tu);
            H_Norm = dH'*dH;
            % Calculate the cost function
            J(i) = tf*(((x1')*x1 + (x2')*x2)/length(t) + ...
                    0.5*R*(u*u')/length(tu));
            % if dH/du < epslon, exit
            if H_Norm < eps
                break;
            else
                % adjust control for next iteration
                u_old = u;
                u = descend(dH,t,u_old,tu,step);
                % modify also Tinit
            end;
        end
        x = [x1';x2'];
        u = interp1(tu,u,t);
        else error('Algorithm not implemented!');
        end
    end
    toc
end

%% State equations for iterative methods
function ydot = stateEq(t,x,u,Tu)
    u = interp1(Tu,u,t); % Interpolate the control at time t
    T = x(end);
    A = [zeros(3),eye(3); zeros(3,6)];
    B = [zeros(3);eye(3)];
    x_dot = T*(A*x(1:6) + B*u);
    ydot = [x_dot; 0];
end

% Costate equations
function pdot = costateEq(t,p,u,tu,x,tx)
    T = p(end);
    x = interp1(tx,x,t); % Interpolate the state variables
    u = interp1(tu,u,t); % Interpolate the control
    A = [zeros(3),eye(3); zeros(3,6)];
    dp = -A' * p(1:6);
    pdot = T*[dp; 0];
end

% Partial derivative of H with respect to u
function dH = Hu(x,p,tx,u,tu)
    % interpolate the control
    u = interp1(tu,u,tx);
    R = 1;
    B = [zeros(3);eye(3)];
    dH = R*u + B'*p;
end

% Adjust the control
function u_new = descend(pH,tx,u,tu,step)
    % interpolate dH/du
    pH = interp1(tx,pH,tu);
    u_new = u - step*pH;
end

%% BVP functions (indirect method)

% -------------------------------------------------------------------------
% ODE's of augmented states
%
function ydot = ode(t,y)
    
    T = y(end);
    x = y(1:6);
    p = y(7:12);
    A = [zeros(3),eye(3); zeros(3,6)];
    B = [zeros(3);eye(3)];
    R = 1;
    u = -R \ (B' * p);
    x_dot = A*x + B*u;
    p_dot = -A' * p;
    ydot = T*[x_dot; p_dot; 0];
end

% -------------------------------------------------------------------------
% boundary conditions: 
%
function res = bc(y0,yf,ballTime,ballPred,robotInit,racketVel)
    % ball at time T
    posRobotInit = robotInit(1:3);
    velRobotInit = robotInit(4:6);
    T = yf(end);
    g = -9.8;
    % Hamiltonian at time T
    p = yf(7:12);
    x = yf(1:6);
    A = [zeros(3),eye(3); zeros(3,6)];
    B = [zeros(3);eye(3)];
    R = 1;
    HatT = -0.5*p'*B*(R\(B'*p)) + p'*A*x;
    ballStateAtT = interp1(ballTime,ballPred',T)';
    ballPosAtT = ballStateAtT(1:3);
    racketVelAtT = interp1(ballTime,racketVel',T)';
    racketAcc = diff(racketVel')./repmat(diff(ballTime'),1,size(racketVel,1));
    racketAccAtT = interp1(ballTime(1:end-1),racketAcc,T)';
    ballVelAtT = ballStateAtT(4:6);
    %bT = posBallInit + velBallInit*T + 0.5*T^2*[0;0;g];
    %dbdT = velBallInit + T*[0;0;g];
    % 6 initial conditions
    res = [ y0(1:3) - posRobotInit;
            y0(4:6) - velRobotInit;
            yf(1:3) - ballPosAtT; 
            yf(4:6) - racketVelAtT;
            %yf(10:12) - 0;
            HatT - yf(7:9)' * ballVelAtT - yf(10:12)' * racketAccAtT];
end
