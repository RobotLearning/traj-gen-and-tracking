%% Three link planar revolute arm inverse kinematics for each joint

% Example taken from Siciliano & Sciavicco, pg. 92-93
% Orientation has to be supplied also as an angle phi (vector)

function q = RRRInvKinematics(x,phi,PAR)

% x is the cartesian coordinate of the endeffector
l1 = PAR.link1.length;
l2 = PAR.link2.length;
l3 = PAR.link3.length;

px = x(1,:) - l3*cos(phi);
py = x(2,:) - l3*sin(phi);

D = (px.^2 + py.^2 - l1^2 - l2^2) / (2*l1*l2);
theta2 = atan2(sqrt(1 - D.^2),D);
theta1 = atan2(py,px) - atan2(l2*sin(theta2),(l1+l2*cos(theta2)));
theta3 = phi - theta2 - theta1;

q = [theta1; theta2; theta3];




