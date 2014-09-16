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
    k1 = h * fun(i,x(:,i),u_trj(:,i),PAR,false);
    % do linear interpolation to get betwixt-u
    x_k1 = x(:,i) + k1/2;
    k2 = h * fun(i+1/2,x_k1,u_trj(:,i),PAR,false);
    x_k2 = x(:,i) + k2/2;
    k3 = h * fun(i+1/2,x_k2,u_trj(:,i),PAR,false);
    x_k3 = x(:,i) + k3;
    k4 = h * fun(i+1,x_k3,u_trj(:,i),PAR,false);
    x(:,i+1) = x(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
