%% Two link planar revolute arm kinematics for each joint

% Example taken from Robot Modelling and Control book, 2006, Spong,
% Hutchingson, Vidyasagar
% pg. 260-262 and 290

function [x1,x2] = RRplanarKinematics(q,PAR)

theta1 = q(1,:);
theta2 = q(2,:);
l1 = PAR.link1.length;
l2 = PAR.link2.length;

x1 = [l1 * cos(theta1); l1 * sin(theta1)];
x2 = [l1 * cos(theta1) + l2 * cos(theta1 + theta2);
      l1 * sin(theta1) + l2 * sin(theta1 + theta2)];


