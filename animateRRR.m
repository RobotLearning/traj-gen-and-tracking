%% creates an animation for the three-link RRR planar arm

function animateRRR(x1,x2,x3,X_des)

% position of the base
%shift = [0; 0.9];
shift = [0;0];

% rotate everything by 90 degrees
R = [0 1; -1 0];
x1 = R * x1;
x2 = R * x2;
x3 = R * x3;
X_des = R * X_des;

x1 = x1 + repmat(shift,1,length(x1));
x2 = x2 + repmat(shift,1,length(x2));
x3 = x3 + repmat(shift,1,length(x3));

% extract desired trajectory
x_des = X_des(1,:) + shift(1);
y_des = X_des(2,:) + shift(2);

% x1 and x2 are the cartesian coordinates of joint 2 and end-effector
link1_x = [shift(1) x1(1,1)];
link1_y = [shift(2) x1(2,1)];
link2_x = [x1(1,1),x2(1,1)];
link2_y = [x1(2,1),x2(2,1)];
link3_x = [x2(1,1),x3(1,1)];
link3_y = [x2(2,1),x3(2,1)];
hf = figure('color','white');
axis equal; axis auto;

h1 = line(link1_x, link1_y, 'color', [.4 .4 .4],'LineWidth',4);
h2 = line(link2_x, link2_y, 'color', [.4 .4 .4],'LineWidth',4);
h3 = line(link3_x, link3_y, 'color', [.4 .4 .4],'LineWidth',4);
title('Two-link planar arm');
hold on;
h4 = scatter(shift(1),shift(2),100,'r','LineWidth',4);
h5 = scatter(x1(1,1),x1(2,1),50,'b','LineWidth',4);
h6 = scatter(x2(1,1),x2(2,1),50,'b','LineWidth',4);
h7 = scatter(x3(1,1),x3(2,1),50,'b','LineWidth',4);
%h6 = plot(x_des,y_des,'-.','color', [.6 .6 .6], 'LineWidth',1);
h8 = scatter(x_des(1),y_des(1),10,'r','LineWidth',2); % for ball tracking

for j = 2:length(x1)
     link1_x(2) = x1(1,j);
     link1_y(2) = x1(2,j);
     link2_x(1) = x1(1,j);
     link2_y(1) = x1(2,j);
     link2_x(2) = x2(1,j);
     link2_y(2) = x2(2,j);
     link3_x(1) = x2(1,j);
     link3_y(1) = x2(2,j);
     link3_x(2) = x3(1,j);
     link3_y(2) = x3(2,j);
     
     set(h1,'XData',link1_x)
     set(h1,'YData',link1_y)
     set(h2,'XData',link2_x)
     set(h2,'YData',link2_y)
     set(h3,'XData',link3_x)
     set(h3,'YData',link3_y)
     set(h5,'XData',link1_x(2));
     set(h5,'YData',link1_y(2));
     set(h6,'XData',link2_x(2));
     set(h6,'YData',link2_y(2));
     set(h7,'XData',link3_x(2));
     set(h7,'YData',link3_y(2));
     
     % for ball tracking
     set(h8,'XData',x_des(j));
     set(h8,'YData',y_des(j));
     
     % connect the resulting end-effector trajectory
     traj_x = [x3(1,j-1),x3(1,j)];               
     traj_y = [x3(2,j-1),x3(2,j)];              
     h9 = line(traj_x, traj_y, 'color', [.4 .4 .5],'LineWidth',1);
     
     drawnow;
     pause(0.005);
end

%hold off;