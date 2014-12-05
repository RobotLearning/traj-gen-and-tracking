%% Two link planar revolute arm inverse kinematics for each joint

% Example taken from Robot Modelling and Control book, 2006, Spong,
% Hutchingson, Vidyasagar
% pg. 260-262 and 290

function q = RRInvKinematics(x,PAR)

% x is the cartesian coordinate of the endeffector
l1 = PAR.link1.length;
l2 = PAR.link2.length;

D = (x(1,:).^2 + x(2,:).^2 - l1^2 - l2^2) / (2*l1*l2);
theta2 = atan(sqrt(1 - D.^2)./D);
theta1 = atan(x(2,:)./x(1,:)) - atan(l2*sin(theta2)./(l1+l2*cos(theta2)));

q = [theta1; theta2];




