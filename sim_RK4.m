function x = sim_RK4(t,x,u_trj,CON,PAR,fun)

% simulates the trajectory of all states using classical Runge-Kutta (RK4)
% the dynamics used is that given by the function handle fun

N = length(t) - 1;
h = t(2) - t(1);
% extract constraints
x_cnstr = CON.state.x.max;
y_cnstr = CON.state.y.max;
w1_cnstr = CON.wheel1.u.max;
w2_cnstr = CON.wheel2.u.max;

for i = 1:N
    % get trajectory of states
    % using classical Runge-Kutta method (RK4)
    k1 = h * fun(i,x(:,i),PAR,u_trj(:,i),false);
    % do linear interpolation to get betwixt-u
    x_k1 = x(:,i) + k1/2;
    u_k_1 = interp1(t, u_trj(1,:), (t(i)+t(i+1))/2, 'linear');
    u_k_2 = interp1(t, u_trj(2,:), (t(i)+t(i+1))/2, 'linear');
    u_k = [u_k_1; u_k_2];
    k2 = h * fun(i+1/2,x_k1,PAR,u_k,false);
    x_k2 = x(:,i) + k2/2;
    k3 = h * fun(i+1/2,x_k2,PAR,u_k,false);
    x_k3 = x(:,i) + k3;
    k4 = h * fun(i+1,x_k3,PAR,u_trj(:,i+1),false);
    x(:,i+1) = x(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    % contraint checking
    if (abs(x(1,i+1)) > x_cnstr || abs(x(2,i+1)) > y_cnstr || ...
        abs(u_trj(1,i+1)) > w1_cnstr || abs(u_trj(2,i+1)) > w2_cnstr)
        error('Constraints violated. Aborting...');
    end
end
