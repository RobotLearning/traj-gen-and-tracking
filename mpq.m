%% Solving Minimum principle for robot kinematics
% Transversality conditions are applied

function [t,x,u,J] = mpq(robotInit,ballTime,ballPred,desVel,jac,kin,djac)

    Tinit = 0.5;
    R = 1;
    momentumInit = zeros(4,1);
    tic;
    % Initial guess for the solution
    solinit = bvpinit(linspace(0,1),[robotInit;momentumInit;Tinit]);
    options = bvpset('Stats','on','RelTol',1e-1);
    bcfull = @(y0,yf) bc(y0,yf,ballTime,ballPred,robotInit,desVel,jac,kin,djac);        
    sol = bvp4c(@ode, bcfull, solinit, options);
    y = sol.y;
    t = y(end)*sol.x;
    y = sol.y;
    x = y(1:4,:);
    % Calculate u from lambda
    p = y(5:8,:);
    B = [zeros(3);eye(3)];
    u = -R \ (B' * p);
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
    
    T = y(end);
    x = y(1:4);
    p = y(5:8);
    A = [zeros(2),eye(2); zeros(2,4)];
    B = [zeros(2);eye(2)];
    R = 1;
    u = -R \ (B' * p);
    x_dot = A*x + B*u;
    p_dot = -A' * p;
    ydot = T*[x_dot; p_dot; 0];
end

% -------------------------------------------------------------------------
% boundary conditions: 
%
function res = bc(y0,yf,ballTime,ballPred,robotInit,racketVel,jac,kin,djac)
    % ball at time T
    posRobotInit = robotInit(1:2);
    velRobotInit = robotInit(3:4);
    T = yf(end);
    g = 9.8;
    % Hamiltonian at time T
    p = yf(5:8);
    q = yf(1:4);
    A = [zeros(2),eye(2); zeros(2,4)];
    B = [zeros(2);eye(2)];
    R = 1;
    HatT = -0.5*p'*B*(R\(B'*p)) + p'*A*q;
    ballStateAtT = interp1(ballTime,ballPred',T)';
    ballPosAtT = ballStateAtT(1:2);
    ballAccAtT = [g;0];
    racketVelAtT = interp1(ballTime,racketVel',T)';
    %racketAcc = diff(racketVel')./repmat(diff(ballTime'),1,size(racketVel,1));
    %racketAccAtT = interp1(ballTime(1:end-1),racketAcc,T)';
    ballVelAtT = ballStateAtT(3:4);
    %bT = posBallInit + velBallInit*T + 0.5*T^2*[0;0;g];
    %dbdT = velBallInit + T*[0;0;g];
    [~,xAtT] = kin(yf(1:2));
    xdAtT = jac(yf(1:2))*yf(3:4);
    % form the A matrix
    Am = [-ballVelAtT', ballAccAtT';
         jac(q(1:2))', djac(q(1:2),q(3:4))';
         zeros(2,2),jac(q(1:2))];
    pars = [HatT;p];
    M = [Am,pars]; % matrix for ensuring transversality
    % 4 initial conditions + 4 terminal cond + 1 time cond from transv.
    res = [ y0(1:2) - posRobotInit;
            y0(3:4) - velRobotInit;
            xAtT - ballPosAtT; 
            xdAtT - racketVelAtT;
            det(M)];
end
