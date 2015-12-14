%% Here we state the boundary conditions to find p(T),x(T),T
% Can be used to solve the bvp with fsolve or lsqnonlin in MATLAB
% The last condition is the transversality condition to find T

function F = bv1(x,PAR)

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

[~,xc] = kin(q); %xc = cart coord.
J = jac(q);
b0 = [b1;b2];
v0 = [v1;v2];
vT = [v1 + g*T;v2];
acc = [g;0];
%vdes = -v;

F = [xc - (b0 + v0*T + (1/2)*acc*T^2);
     nu - 0;
     %J * qdot - vdes;
     %H - coeff * (p1*l2*cos(q1+q2)*v1 + p1*l2*sin(q1+q2)*v2 - ...
     %     p2*l1*v1*cos(q1) - p2*l2*v1*cos(q1+q2) - p2*l1*v2*sin(q1) - p2*l2*v2*sin(q1+q2))];
     H - mu'*(J\vT)];
%{
F = [l1*cos(q1) + l2*cos(q1+q2)- b1 - v1*T;
     l1*sin(q1) + l2*sin(q1+q2) - b2 - v2*T;
     (-l1*sin(q1)-l2*sin(q1+q2))*qdot1 - l2*sin(q1+q2)*qdot2 + v1;
     (-l1*cos(q1)-l2*cos(q1+q2))*qdot1 - l2*cos(q1+q2)*qdot2 + v2;
     H - (denum * p1*l2*cos(q1+q2)*v1 + p1*l2*sin(q1+q2)*v2 - ...
          p2*l1*v1*cos(q1) - p2*l2*v1*cos(q1+q2) - p2*l1*v2*sin(q1) - p2*l2*v2*sin(q1+q2))];
%}     


% b1 = par(1);
% b2 = par(2);
% b3 = par(3);
% v1 = par(4);
% v2 = par(5);
% v3 = par(6);
% g = -9.8;
% 
% p1 = x(1);
% p2 = x(2);
% p3 = x(3);
% p4 = x(4);
% p5 = x(5);
% p6 = x(6);
% T = x(7);
% 
% F = [1/6*p1*T^3 + 1/2*p4*T^2 - v1*T - b1; ... % desired ball pos at T
%         1/6*p2*T^3 + 1/2*p5*T^2 - v2*T - b2; ...
%         1/6*p3*T^3 + 1/2*p6*T^2 - v3*T - 1/2*g*T^2 - b3; ...
%         1/2*p1*T^2 + p4*T + v1; ... % desired racket velocities at T
%         1/2*p2*T^2 + p5*T + v2; ...
%         1/2*p3*T^2 + p6*T + v3 + g*T; ...
%          p1*v1 + p2*v2 + p3*(v3+g*T) + (p6+p3*T)*g + ...
%          +1/2*(p4^2 + p5^2 + p6^2) ]; 


%F = [-p1*T + x1 - b1 - v1*T; ...
%     -p2*T + x2 - b2 - v2*T; ...
%      1/2*(p1^2 + p2^2) + p1*v1 + p2*v2];