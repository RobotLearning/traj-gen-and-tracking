% Dynamic motor primitives

function discreteAccelerationDMP

g = 2;

% alpha_z and beta_z are time constants
alpha_z = 1;
beta_z = 0.25;
alpha_v = 1;
beta_v = 0.25;
alpha_g = 0.5;

% temporal scaling factor
tau = 0.5;

% start from z,y nonzero
z0 = 0;
y0 = 0;
r0 = 0;
v0 = 0;
x0 = 0;
X0 = [z0;y0;r0;v0;x0];

% forcing function parameters
f.w = [0.5, 0.5, 0.5];
f.c = [2, 1, 0];
f.h = [1, 1, 1];

% time evolution
tspan = [0 10];

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4]);
PAR.time_const = [alpha_z, beta_z, alpha_v, beta_v, alpha_g];
PAR.scale = tau;
PAR.goal = g;
PAR.forcing = f;
[t,X] = ode45(@(t,X)dmp1(t,X,PAR),tspan,X0,options);
z = X(:,1);
y = X(:,2);
r = X(:,3);
v = X(:,4);
x = X(:,5);
plot(t,z,'-',t,y,'-.',t,v,t,x,t,g*ones(1,length(t)),'r.');
legend('state velocity z = ydot','state y','phase velocity',...
       'phase','goal state');

end

function dX = dmp1(t,X,PAR)

alpha_z = PAR.time_const(1);
beta_z = PAR.time_const(2);
alpha_v = PAR.time_const(3);
beta_v = PAR.time_const(4);
alpha_g = PAR.time_const(5);
tau = PAR.scale;
g = PAR.goal;
fpar = PAR.forcing;

A = [-alpha_z/tau, -alpha_z*beta_z/tau, alpha_z*beta_z/tau, zeros(1,2);
     1/tau, 0, 0, zeros(1,2);
     zeros(1,2), -alpha_g/tau, zeros(1,2);
     zeros(1,3), -alpha_v/tau, -alpha_v*beta_v/tau;
     zeros(1,3), 1/tau, 0];

f = forcing(X(5),X(4),g,fpar.w,fpar.h,fpar.c);
% forcing function acts on the accelerations
B = [f/tau; 0; alpha_g * g/tau; alpha_v*beta_v*g/tau; 0];

dX = A*X + B;

end

% basis functions are unscaled gaussians
function phi = basis(x,g,h,c)
phi = exp(-h * (x/g - c)^2);
end

% forcing function to drive nonlinear system dynamics
function f = forcing(x,v,g,w,h,c)

N = length(w);
f = 0;
scale = 0;
for i = 1:N
    f = f + basis(x,g,h(i),c(i))*w(i)*v;
    scale = scale + basis(x,g,h(i),c(i));
end

f = f/scale;
end