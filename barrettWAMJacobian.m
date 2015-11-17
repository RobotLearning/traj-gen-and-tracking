% Barrett WAM Jacobian
% used to show racket velocities in cartesian space
%
% Function taken from SL: 
% shared/barrett/include/SL_kinematics_body.h
%
% Called from the kinematics method of Barrett WAM class
%
%
function xd = barrettWamJacobian(Q,PAR)

NDOF = 7;
% system states are X   =   (q(1),...q(7),qd(1),...,qd(7));
q = Q(1:NDOF);
qd = Q(NDOF+1:2*NDOF);

% definitions
ZSFE = 0.346;              %!< z height of SAA axis above ground
ZHR = 0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB = 0.045;              %!< elbow y offset
ZEB = 0.045;              %!< elbow z offset
YWR = -0.045;              %!< elbow y offset (back to forewarm)
ZWR = 0.045;              %!< elbow z offset (back to forearm)
ZWFE = 0.255;              %!< forearm length (minus 4.5cm)

% extract parameters
link0 = PAR.link0;
links = PAR.links;
eff = PAR.eff;
uex = PAR.uex;
uex0 = PAR.uex0;
basec = PAR.basec;
baseo = PAR.baseo;