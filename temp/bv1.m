function F = bv1(x,PAR)

% Solve symbolically to compare
p = x(1:end-1);
dim = length(p)/2;
% divide as mu and nu
mu = p(1:dim);
nu = p(dim+1:end);
T = x(end);

l1 = PAR.l1;
l2 = PAR.l2;
b1 = PAR.b1;
b2 = PAR.b2;
v1 = PAR.v1;
v2 = PAR.v2;
q10 = PAR.q1;
q20 = PAR.q2;
kin = PAR.kin; % kinematics function
jac = PAR.jac; % jacobian

q = -(1/6)*mu*T^3 - (1/2)*nu*T^2 + q10;
qdot = (-1/2)*mu*T^2 - nu*T;

H = -1/2*(nu'*nu);
%coeff = 1/(l1*l2*(cos(q1)*sin(q1+q2) - sin(q1)*cos(q1+q2)));

[~,xc] = kin(q); %xc = cart coord.
J = jac(q);
b = [b1;b2];
v = [v1;v2];
vdes = -v;

F = [xc - b - v*T;
     J * qdot - vdes;
     %H - coeff * (p1*l2*cos(q1+q2)*v1 + p1*l2*sin(q1+q2)*v2 - ...
     %     p2*l1*v1*cos(q1) - p2*l2*v1*cos(q1+q2) - p2*l1*v2*sin(q1) - p2*l2*v2*sin(q1+q2))];
     H - mu'*(J\v)];
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