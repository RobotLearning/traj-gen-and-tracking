% Dynamic motor primitives

function discreteAccelerationDMP


% time constant of the canonical system
alpha_x = 1;
% temporal scaling factor
tau = 1;
% start from nonzero phase
x0 = 1;
% time evolution
tspan = [0 10];

% goal state
DMP(1).g = 1;
DMP(2).g = 5;
% alpha_g and beta_g are time constants
DMP(1).alpha_g = 1;
DMP(1).beta_g = 0.5;
DMP(1).alpha_x = alpha_x;
DMP(2).alpha_g = 1;
DMP(2).beta_g = 0.5;
DMP(2).alpha_x = alpha_x;

% start from z,y nonzero
z0 = 0;
y0 = 0;
DMP(1).X0 = [z0;y0;x0];
z0 = 0;
y0 = 1;
DMP(2).X0 = [z0;y0;x0];

% forcing function parameters
DMP(1).forcing.w = [1, .1, .5];
DMP(1).forcing.c = [2, 1, 0];
DMP(1).forcing.h = [1, .1, 1];
DMP(2).forcing.w = [1, 1, 1];
DMP(2).forcing.c = [2, 1, 3];
DMP(2).forcing.h = [1, .1, 3];


options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
DMP(1).scale = tau;
DMP(2).scale = tau;
[t1,DMP(1).X] = ode45(@(t,X)funDiscrete(t,X,DMP(1)),tspan,DMP(1).X0,options);
[t2,DMP(2).X] = ode45(@(t,X)funDiscrete(t,X,DMP(2)),tspan,DMP(2).X0,options);

% plotting each dmp trajectory
z = DMP(1).X(:,1);
y = DMP(1).X(:,2);
x = DMP(1).X(:,3);
figure(1);
plot(t1,z,'-',t1,y,'-.',t1,x,t1,DMP(1).g*ones(1,length(t1)),'r.');
legend('state velocity ydot','state y','phase','goal state');
z = DMP(2).X(:,1);
y = DMP(2).X(:,2);
x = DMP(2).X(:,3);
figure(2);
plot(t2,z,'-',t2,y,'-.',t2,x,t2,DMP(2).g*ones(1,length(t2)),'r.');
legend('state velocity ydot','state y','phase','goal state');

end




