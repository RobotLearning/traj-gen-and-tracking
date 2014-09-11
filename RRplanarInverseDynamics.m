% Two-link planar revolute arm inverse dynamics for evolving system
% trajectories

% Example taken from Robot Modelling and Control book, 2006, Spong,
% Hutchingson, Vidyasagar
% pg. 260-262 and 290

function Qd = RRplanarInverseDynamics(Q,u,PAR)

% system states are X = [q(1),q(2),qd(1),qd(2)];
q = Q(1:2);
qd = Q(3:4);

g = PAR.const.g;
m1 = PAR.link1.mass;
m2 = PAR.link2.mass;
l1 = PAR.link1.length;
l2 = PAR.link2.length;
l_c1 = PAR.link1.centre.dist;
l_c2 = PAR.link2.centre.dist;
I1 = PAR.link1.inertia;
I2 = PAR.link2.inertia;
J_m1 = PAR.link1.motor.inertia;
J_m2 = PAR.link2.motor.inertia;
r_1 = PAR.link1.motor.gear_ratio;
r_2 = PAR.link2.motor.gear_ratio;

% inertia matrix due to the rotors
J = [r_1^2 * J_m1, 0;
     0, r_2^2 * J_m2];

jac1 = [-l_c1 * sin(q(1)), 0;
        l_c1 * cos(q(1)), 0;
        0, 0];
jac2 = [-l1*sin(q(1)) - l_c2*sin(q(1)+q(2)), -l_c2*sin(q(1)+q(2));
        l1*cos(q(1)) + l_c2*cos(q(1)+q(2)), l_c2*cos(q(1)+q(2));
        0, 0];
% Translational part of the kinetic energy
T = m1 * (jac1' * jac1) + m2 * (jac2' * jac2);
% Angular part of the kinetic energy
R = I1 * [1,0;0,0] + I2 * ones(2);

% The inertia matrix
D = T + R;
% d11 = m1*l_c1^2 + m2*(l1^2 + l_c2^2 + 2*l1*l_c2*cos(q(2))) + I1 + I2;
% d12 = m2*(l_c2^2 + l1*l_c2*cos(q(2))) + I2;
% d22 = m2*l_c2^2 + I2;
% D = [d11, d12;
%      d12, d22];
M = D + J;

% Coriolis and centrifugal generalized forces
c111 = 0;
c121 = -m2*l1*l_c2*sin(q(2));
c211 = c121;
c221 = -m2*l1*l_c2*sin(q(2));
c112 = m2*l1*l_c2*sin(q(2));
c122 = 0;
c212 = 0;
c222 = 0;

Cpr = [c111, c211;
       c121, c221;
       c112, c212;
       c122, c222] * qd;
C = [Cpr(1), Cpr(2); Cpr(3), Cpr(4)];

% friction matrix is assumed to be zero
B = 0;

% potential related (conservative) forces
g1 = (m1*l_c1 + m2*l1)*g*cos(q(1)) + m2*l_c2*g*cos(q(1)+q(2));
g2 = m2*l_c2*g*cos(q(1)+q(2));
G = [g1; g2];

% control inputs are the scaled voltages sent to the rotors
% u_k = r_k * K_mk / R_k * V_k

% construct state-dependent drift term
Ax = [zeros(2), eye(2);
      zeros(2), -C];

% construct control matrix
Bx = [zeros(2); eye(2)];

% state-independent drift term
Cx = [zeros(2,1); -G];

Mbig = [eye(2), zeros(2); zeros(2), M];
Qd = Mbig \ (Ax*Q + Bx*u + Cx);
