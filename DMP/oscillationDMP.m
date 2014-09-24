% Dynamic motor primitives
% TODO: doesn't track the goal state!

function oscillationDMP

% set point of the goal
g = 1;
% amplitude of the goal
A = 1;

% alpha_z and beta_z are time constants
alpha_z = 1;
beta_z = 1;

% temporal scaling factor
tau = 1;

% start from z,y nonzero
z0 = 0;
y0 = 0;
phi0 = 0;
X0 = [z0;y0;phi0];

% forcing function parameters
f.w = [1, .1, .5];
f.c = [2, 1, 5];
f.h = [1, 1, 1];

% time evolution
tspan = [0 20];

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
PAR.time_const = [alpha_z, beta_z];
PAR.scale = tau;
PAR.goal = [g, A];
PAR.forcing = f;
[t,X] = ode45(@(t,X)funRhythmic(t,X,PAR),tspan,X0,options);
z = X(:,1);
y = X(:,2);
phi = X(:,3);
plot(t,z,'-',t,y,'-.',t,sin(phi),t,g*ones(1,length(t)),'r.');
legend('state velocity z = ydot','state y',...
       'phase','goal state');

end