%% Solving Minimum principle 

function [Tx,X,Tu,u,J] = mp()

    t0 = 0; tf = 0.78;
    N = 100;
    Tu = linspace(t0,tf,N);
    X(1:2,1) = [0.05;0];
    max_iteration = 100;
    u = zeros(1,N);
    
    % Initial guess for the solution
    solinit = bvpinit(linspace(t0,tf,N), [0 0 0.5 0.5]);
    options = bvpset('Stats','on','RelTol',1e-1);
    global R;
    R = 0.1;
    sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
    t = sol.x;
    y = sol.y;
    % Calculate u(t) from x1,x2,p1,p2
    ut = (y(3,:).*(y(1,:) + 1/4))/(2*0.1);
    n = length(t);
    % Calculate the cost
    J = tf*(y(1,:)*y(1,:)' + y(2,:)*y(2,:)' + ...
    ut*ut'*0.1)/n;
    % arrange output
    Tx = t;
    X = y;
    u = ut;
    
    
    %% Steepest descent
%     J = zeros(1,max_iteration);
%     eps = 1e-3;
%     step = 0.2;
%     options = [];
%     for i = 1:max_iteration
%         % 1) start with assumed control u and move forward
%         initx = X(1:2,1);
%         [Tx,X] = ode45(@(t,x) stateEq(t,x,u,Tu), [t0 tf], ...
%                        initx, options);
%         % 2) Move backward to get the trajectory of costates
%         x1 = X(:,1); x2 = X(:,2);
%         % calculate initp
%         initp = [0;0];
%         [Tp,P] = ode45(@(t,p) costateEq(t,p,u,Tu,x1,x2,Tx), ...
%         [tf t0], initp, options);
%         p1 = P(:,1);
%         % Important: costate is stored in reverse order. The dimension of
%         % costates may also different from dimension of states
%         % Use interploate to make sure x and p is aligned along the time axis
%         p1 = interp1(Tp,p1,Tx);
%         % Calculate deltaH with x1(t), x2(t), p1(t), p2(t)
%         dH = Hu(x1,p1,Tx,u,Tu);
%         H_Norm = dH'*dH;
%         % Calculate the cost function
%         J(i) = tf*(((x1')*x1 + (x2')*x2)/length(Tx) + ...
%                 0.1*(u*u')/length(Tu));
%         % if dH/du < epslon, exit
%         if H_Norm < eps
%             % Display final cost
%             J(i)
%             break;
%         else
%             % adjust control for next iteration
%             u_old = u;
%             u = descend(dH,Tx,u_old,Tu,step);
%         end;
%     end
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
function dydt = BVP_ode(t,y)
    global R;
    t1 = y(1)+.25;
    t2 = y(2)+.5;
    t3 = exp(25*y(1)/(y(2)+2));
    t4 = 50/(y(1)+2)^2;
    u = y(3)*t1/(2*R);
    dydt = [-2*t1+t2*t3-t2*u
    0.5-y(2)-t2*t3
    -2*y(1)+2*y(3)-y(3)*t2*t4*t3+y(3)*u+y(4)*t2*t4*t3
    -2*y(2)-y(3)*t3+y(4)*(1+t3)];
end

% -----------------------------------------------
% The boundary conditions:
% x1(0) = 0.05, x2(0) = 0, tf = 0.78, p1(tf) = 0, p2(tf) = 0;
%
function res = BVP_bc(ya,yb)
    res = [ ya(1) - 0.05
    ya(2) - 0
    yb(3) - 0
    yb(4) - 0 ];
end
