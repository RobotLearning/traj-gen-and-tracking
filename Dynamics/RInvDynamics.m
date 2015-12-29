% Pendulum Inverse Dynamics
%

function u = RInvDynamics(q,qd,qdd,PAR)

g = PAR.const.g;
m = PAR.link.mass;
l = PAR.link.length;
l_c = PAR.link.centre.dist;
I = PAR.link.inertia;
J_m = PAR.link.motor.inertia;
r = PAR.link.motor.gear_ratio;

% inertia matrix due to the rotors
J = r^2 * J_m;

jac = [-l_c * sin(q(1)); l_c * cos(q(1)); 0];
% Translational part of the kinetic energy
T = m * (jac' * jac);
% Angular part of the kinetic energy
R = I;

% The inertia matrix
D = T + R;
M = D + J;

% Coriolis and centrifugal generalized forces
C = 0;

% friction matrix is assumed to be zero
B = 0;

% potential related (conservative) forces
G = m*l_c*g*cos(q(1));

% control inputs are the scaled voltages sent to the rotors
% u_k = r_k * K_mk / R_k * V_k
u = M*qdd + C*qd + B*qd + G;