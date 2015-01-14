% Pendulum Dynamics
%
% If flag is set to true, return the jacobians of the dynamics f

function [Qd,varargout] = RDynamics(Q,u,PAR,flag)

% system states are X = [q(1),qd(1)];
q = Q(1);
qd = Q(2);

g = PAR.const.g;
m = PAR.link.mass;
l = PAR.link.length;
l_c = PAR.link.centre.dist;
I = PAR.link.inertia;
J_m = PAR.link.motor.inertia;
r = PAR.link.motor.gear_ratio;

% inertia matrix due to the rotor
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

% construct state-dependent drift term
Ax = [0, 1; 0, -C];

% construct control matrix
Bx = [0; 1];

% state-independent drift term
Cx = [0; -G];

Mbig = [1, 0; 0, M];
Qd = Mbig \ (Ax*Q + Bx*u + Cx);

if flag
    
    % Take four differences with h small
    h = 1e-6;
    Qh = repmat(Q,1,2) + h * eye(2);
    Qdh = zeros(2);
    for i = 1:2
        Qdh(:,i) = RDynamics(Qh(:,i),u,PAR,false);
    end
    der = (Qdh - repmat(Qd,1,2)) / h;
    dfdx = [0, 1; der(2,:)];
    dfdu = Mbig \ Bx;
    varargout{1} = dfdx;
    varargout{2} = dfdu;
    
end