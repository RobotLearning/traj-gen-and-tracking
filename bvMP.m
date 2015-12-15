%% Here we state the boundary conditions to find p(T),x(T),T
% Can be used to solve the maximum principle bvp 
% with fsolve or lsqnonlin in MATLAB
% The last condition is the transversality condition to find T

function F = bvMP(x,PAR)

% Solve symbolically to compare
p = x(1:end-1);
dim = length(p)/2;
% divide as mu and nu
mu = p(1:dim);
nu = p(dim+1:end);
T = x(end);

b1 = PAR.b1;
b2 = PAR.b2;
v1 = PAR.v1;
v2 = PAR.v2;
q0 = PAR.q0;
kin = PAR.kin; % kinematics function
jac = PAR.jac; % jacobian
g = PAR.g;

% since we're minimizing accelerations 3rd degree polynomials
mu0 = mu;
nu0 = -nu - mu0*T;
q = (1/6)*mu0*T^3 + (1/2)*nu0*T^2 + q0;
qdot = (1/2)*mu0*T^2 + nu0*T;

H = (1/2)*(mu0*T - nu0)'*(mu0*T + nu0); %0; %-1/2*(nu'*nu);
%coeff = 1/(l1*l2*(cos(q1)*sin(q1+q2) - sin(q1)*cos(q1+q2)));

[~,~,xc] = kin(q); %xc = cart coord.
J = jac(q);
b0 = [b1;b2];
v0 = [v1;v2];
vT = [v1 + g*T;v2];
acc = [g;0];
%vdes = -v;
Dphi = [-vT, J];
r = (eye(dim+1) - Dphi'*pinv(Dphi)')*[H;-mu];

F = [xc - (b0 + v0*T + (1/2)*acc*T^2);
     nu - 0;
     %J * qdot - vdes;
     %H - coeff * (p1*l2*cos(q1+q2)*v1 + p1*l2*sin(q1+q2)*v2 - ...
     %     p2*l1*v1*cos(q1) - p2*l2*v1*cos(q1+q2) - p2*l1*v2*sin(q1) - p2*l2*v2*sin(q1+q2))];
     %H - mu'*(J\vT)];
     r];
