%% Solving Minimum principle for 2d kinematics case
% Transversality conditions are applied

function [t,x,u,J] = mp0(robotInit,ballTime,ballPred)

    Tinit = 0.5;
    R = 1;
    momentumInit = zeros(2,1);
    tic;
    % Initial guess for the solution
    solinit = bvpinit(linspace(0,1),[robotInit;momentumInit;Tinit]);
    options = bvpset('Stats','on','RelTol',1e-1);
    bcfull = @(y0,yf) bc(y0,yf,ballTime,ballPred,robotInit);        
    sol = bvp4c(@ode, bcfull, solinit, options);
    y = sol.y;
    t = y(end)*sol.x;
    y = sol.y;
    x = y(1:2,:);
    % Calculate u from lambda
    p = y(3:4,:);
    u = -R\p;
    % Calculate cost
    T = t(end);
    N = length(t);
    dt = T/N;
    J = 0.5 * dt * trace(u'*R*u);
    toc
end
%% BVP functions (indirect method)

% -------------------------------------------------------------------------
% ODE's of augmented states
%
function ydot = ode(t,y)
    
    R = 1;
    T = y(end);
    x = y(1:2);
    p = y(3:4);
    u = -R \ p;
    x_dot = u;
    ydot = T*[x_dot; 0; 0; 0];
end

% -------------------------------------------------------------------------
% boundary conditions: 
%
function res = bc(y0,yf,ballTime,ballPred,robotInit)
    % ball at time T
    posRobotInit = robotInit(1:2);
    T = yf(end);
    % Hamiltonian at time T
    p = yf(3:4);
    x = yf(1:2);
    R = 1;
    HatT = -0.5*(p'*p);
    ballStateAtT = interp1(ballTime,ballPred',T)';
    ballPosAtT = ballStateAtT(1:2);
    ballVelAtT = ballStateAtT(3:4);
    % 2 initial conditions
    res = [ y0(1:2) - posRobotInit;
            yf(1:2) - ballPosAtT; 
            %yf(10:12) - 0;
            HatT - p' * ballVelAtT];
end
