% Dynamic motor primitives

function oscillationDMP

% set point of the goal
g = 1;
% amplitude of the goal
A = 10;

% alpha_z and beta_z are time constants
alpha_z = 1;
beta_z = 0.25;

% temporal scaling factor
tau = 1;

% start from z,y nonzero
z0 = 0;
y0 = 0;
phi0 = 0;
X0 = [z0;y0;phi0];

% forcing function parameters
f.w = [0.5, 1, 0.5];
f.c = [2, 0, 5];
f.h = [10, 3, 2];

% time evolution
tspan = [0 10];

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
PAR.time_const = [alpha_z, beta_z];
PAR.scale = tau;
PAR.goal = [g, A];
PAR.forcing = f;
[t,X] = ode45(@(t,X)dmp1(t,X,PAR),tspan,X0,options);
z = X(:,1);
y = X(:,2);
phi = X(:,3);
plot(t,z,'-',t,y,'-.',t,A*cos(phi),t,g*ones(1,length(t)),'r.');
legend('state velocity z = ydot','state y',...
       'phase','goal state');

end

function dX = dmp1(t,X,PAR)

alpha_z = PAR.time_const(1);
beta_z = PAR.time_const(2);
tau = PAR.scale;
g = PAR.goal(1);
Amp = PAR.goal(2);
fpar = PAR.forcing;

A = [-alpha_z/tau, -alpha_z*beta_z/tau, 0;
     1/tau, 0, 0;
     0, 0, 0];

f = forcing(X(3),Amp,fpar.w,fpar.h,fpar.c);
% forcing function acts on the accelerations
B = [alpha_z*beta_z*g/tau + f/tau; 0; 1/tau];

dX = A*X + B;

end

% basis functions are unscaled gaussians
function phi = basis(phi,h,c)
phi = exp(h * (cos(phi - c) - 1));
end

% forcing function to drive nonlinear system dynamics
function f = forcing(phi,A,w,h,c)

N = length(w);
f = 0;
scale = 0;
for i = 1:N
    f = f + basis(phi,h(i),c(i))*w(i);
    scale = scale + basis(phi,h(i),c(i));
end

f = f/scale * A;
end