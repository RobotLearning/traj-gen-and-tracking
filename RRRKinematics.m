%% Three link planar revolute arm kinematics for each joint

% Example taken from Robotics, Modelling Planning and Control
% Siciliano Sciavicco, Springer, pg. 70

function [x1,x2,x3,Ahmats] = RRRKinematics(q,PAR)

theta1 = q(1,:);
theta2 = q(2,:);
theta3 = q(3,:);
l1 = PAR.link1.length;
l2 = PAR.link2.length;
l3 = PAR.link3.length;
c1 = cos(theta1);
s1 = sin(theta1);
c12 = cos(theta1 + theta2);
s12 = sin(theta1 + theta2);
c123 = cos(theta1 + theta2 + theta3);
s123 = sin(theta1 + theta2 + theta3);


x1 = [l1 * c1; l1 * s1];
x2 = [l1 * c1 + l2 * c12;
      l1 * s1 + l2 * s12];
x3 = [l1 * c1 + l2 * c12 + l3 * c123;
      l1 * s1 + l2 * s12 + l3 * s123];

Ahmats = zeros(4,4,3,size(theta1,2));
for i = 1:size(theta1,2)  
    R1 = [c1(i),-s1(i),0;
          s1(i),c1(i),0;
          0,0,1];
    R2 = [c12(i),-s12(i),0;
          s12(i),c12(i),0;
          0,0,1];
    R3 = [c123(i),-s123(i),0;
          s123(i),c123(i),0;
          0,0,1];
    p1 = [l1*c1(i);l1*s1(i);0];
    p2 = [l1*c1(i)+l2*c12(i);l1*s1(i)+l2*s12(i);0];
    p3 = [l1*c1(i)+l2*c12(i)+l3*c123(i);
          l1*s1(i)+l2*s12(i)+l3*s123(i);0];
    Ahmats(:,:,1,i) = [R1,p1;zeros(1,3),1];
    Ahmats(:,:,2,i) = [R2,p2;zeros(1,3),1];
    Ahmats(:,:,3,i) = [R3,p3;zeros(1,3),1];
end



