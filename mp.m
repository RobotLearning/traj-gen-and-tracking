%% Solving Minimum principle 

function [t,x,u,J] = mp(solve_method)

    t0 = 0; tf = 0.78;
    N = 100;
    x(1:2,1) = [0.05;0];
    max_iteration = 100;
    u = zeros(1,N);
    tic;
    if strcmp(solve_method,'BVP')
        tic;
        % Initial guess for the solution
        solinit = bvpinit(linspace(t0,tf,N), [0 0 0.5 0.5]);
        options = bvpset('Stats','on','RelTol',1e-1);
        R = 0.1;
        sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
        t = sol.x;
        y = sol.y;
        x = y(1:2,:);
        p = y(3:4,:);
        % Calculate u(t) from x1,x2,p1,p2
        u = (p(1,:).*(x(1,:) + 1/4))/(2*R);
        % Calculate the cost
        dt = tf/length(t);
        J = dt * (x(1,:)*x(1,:)' + x(2,:)*x(2,:)' + u*R*u');
    else    
        %% Steepest descent
        J = zeros(1,max_iteration);
        eps = 1e-3;
        step = 0.2;
        options = [];
        tu = linspace(t0,tf,N);
        for i = 1:max_iteration
            % 1) start with assumed control u and move forward
            initx = x(1:2,1);
            [t,x] = ode45(@(t,x) stateEq(t,x,u,tu), [t0 tf], ...
                           initx, options);
            % 2) Move backward to get the trajectory of costates
            x1 = x(:,1); x2 = x(:,2);
            % calculate initp
            initp = [0;0];
            [tp,P] = ode45(@(tp,p) costateEq(tp,p,u,tu,x1,x2,t), ...
            [tf t0], initp, options);
            p1 = P(:,1);
            % Important: costate is stored in reverse order. The dimension of
            % costates may also different from dimension of states
            % Use interploate to make sure x and p is aligned along the time axis
            p1 = interp1(tp,p1,t);
            % Calculate deltaH with x1(t), x2(t), p1(t), p2(t)
            dH = Hu(x1,p1,t,u,tu);
            H_Norm = dH'*dH;
            % Calculate the cost function
            J(i) = tf*(((x1')*x1 + (x2')*x2)/length(t) + ...
                    0.1*(u*u')/length(tu));
            % if dH/du < epslon, exit
            if H_Norm < eps
                break;
            else
                % adjust control for next iteration
                u_old = u;
                u = descend(dH,t,u_old,tu,step);
            end;
        end
        x = [x1';x2'];
        u = interp1(tu,u,t);
    end
    elapsedTime = toc;
    fprintf('MP took %f sec \n', elapsedTime);
end

% State equations
function dx = stateEq(t,x,u,Tu)
    dx = zeros(2,1);
    u = interp1(Tu,u,t); % Interpolate the control at time t
    dx(1) = -2*(x(1) + 0.25) + (x(2) + 0.5)*exp(25*x(1)/(x(1)+2)) ...
            -(x(1) + 0.25).*u;
    dx(2) = 0.5 - x(2) -(x(2) + 0.5)*exp(25*x(1)/(x(1)+2));
end

% Costate equations
function dp = costateEq(t,p,u,Tu,x1,x2,xt)
    dp = zeros(2,1);
    x1 = interp1(xt,x1,t); % Interpolate the state variables
    x2 = interp1(xt,x2,t);
    u = interp1(Tu,u,t); % Interpolate the control
    dp(1) = p(1).*(u + exp((25*x1)/(x1 + 2)).*((25*x1)/(x1 + 2)^2 - ...
    25/(x1 + 2))*(x2 + 1/2) + 2) - ...
    2*x1 - p(2).*exp((25*x1)/(x1 + 2))*((25*x1)/(x1 + 2)^2 - ...
    25/(x1 + 2))*(x2 + 1/2);
    dp(2) = p(2).*(exp((25*x1)/(x1 + 2)) + 1) - ...
    p(1).*exp((25*x1)/(x1 + 2)) - 2*x2;
end

% Partial derivative of H with respect to u
function dH = Hu(x1,p1,tx,u,Tu)
    % interpolate the control
    u = interp1(Tu,u,tx);
    R = 0.1;
    dH = 2*R*u - p1.*(x1 + 0.25);
end

% Adjust the control
function u_new = descend(pH,tx,u,tu,step)
    % interpolate dH/du
    pH = interp1(tx,pH,tu);
    u_new = u - step*pH;
end

%% BVP functions

%------------------------------------------------
% ODE’s for states and costates
%
function dxdt = BVP_ode(t,x)
    R = 0.1;
    t1 = x(1) + .25;
    t2 = x(2) + .5;
    t3 = exp(25*x(1)/(x(2) + 2));
    t4 = 50/(x(1) + 2)^2;
    u = x(3)*t1/(2*R);
    dxdt = [-2*t1 + t2*t3 - t2*u;
            0.5-x(2) - t2*t3;
            -2*x(1) + 2*x(3) - x(3)*t2*t4*t3 + x(3)*u + x(4)*t2*t4*t3;
            -2*x(2) - x(3)*t3 + x(4)*(1+t3)];
end

% -----------------------------------------------
% The boundary conditions:
% x1(0) = 0.05, x2(0) = 0, tf = 0.78, p1(tf) = 0, p2(tf) = 0;
%
function res = BVP_bc(xa,xb)
    res = [ xa(1) - 0.05;
            xa(2) - 0;
            xb(3) - 0;
            xb(4) - 0];
end
