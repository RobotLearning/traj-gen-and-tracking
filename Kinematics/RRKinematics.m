%% Two link planar revolute arm kinematics for each joint

% Example taken from Robot Modelling and Control book, 2006, Spong,
% Hutchingson, Vidyasagar
% pg. 260-262 and 290

function [x1,x2,Ahmats] = RRKinematics(q,PAR)

theta1 = q(1,:);
theta2 = q(2,:);
l1 = PAR.link1.length;
l2 = PAR.link2.length;
c1 = cos(theta1);
s1 = sin(theta1);
c12 = cos(theta1 + theta2);
s12 = sin(theta1 + theta2);

x1 = [l1 * c1; l1 * s1];
x2 = [l1 * c1 + l2 * c12;
      l1 * s1 + l2 * s12];

for i = 1:length(theta1)  
    R1 = [c1(i),-s1(i),0;
          s1(i),c1(i),0;
          0,0,1];
    R2 = [c12(i),-s12(i),0;
          s12(i),c12(i),0;
          0,0,1];
    p1 = [l1 * c1(i); l1 * s1(i); 0];
    p2 = [l1 * c1(i) + l2 * c12(i); l1 * s1(i) + l2 * s12(i); 0];
    Ahmats(:,:,1,i) = [R1, p1; zeros(1,3), 1];
    Ahmats(:,:,2,i) = [R2, p2; zeros(1,3), 1];
end
