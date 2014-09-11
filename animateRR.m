%% creates an animation for the two-link RR planar arm

function animateRR(x1,x2,X_des)

% extract desired trajectory
x_des = X_des(1,:);
y_des = X_des(2,:);

% x1 and x2 are the cartesian coordinates of joint 2 and end-effector

link1_x = [0 x1(1,1)];
link1_y = [0 x1(2,1)];
link2_x = [x1(1,1),x2(1,1)];
link2_y = [x1(2,1),x2(2,1)];
hf = figure('color','white');
axis manual, axis equal;
h1 = line(link1_x, link1_y, 'color', [.4 .4 .8],'LineWidth',4);
h2 = line(link2_x, link2_y, 'color', [.4 .4 .8],'LineWidth',4);
title('Two-link planar arm','Color',[.6 0 0]);
hold on;
h3 = scatter(0,0,100,'r','LineWidth',4);
h4 = scatter(x1(1,1),x1(2,1),100,'r','LineWidth',4);
h5 = scatter(x2(1,1),x2(2,1),100,'r','LineWidth',4);
h6 = plot(x_des,y_des,'-r','LineWidth',2);

for j = 2:length(x1)
     link1_x(2) = x1(1,j);
     link1_y(2) = x1(2,j);
     link2_x(1) = x1(1,j);
     link2_y(1) = x1(2,j);
     link2_x(2) = x2(1,j);
     link2_y(2) = x2(2,j);
     
     set(h1,'XData',link1_x)
     set(h1,'YData',link1_y)
     set(h2,'XData',link2_x)
     set(h2,'YData',link2_y)
     set(h4,'XData',link1_x(2));
     set(h4,'YData',link1_y(2));
     set(h5,'XData',link2_x(2));
     set(h5,'YData',link2_y(2));
     
     % connect the resulting end-effector trajectory
     traj_x = [x2(1,j-1),x2(1,j)];               
     traj_y = [x2(2,j-1),x2(2,j)];              
     h7 = line(traj_x, traj_y, 'color', [.4 .4 .5],'LineWidth',2);
     
     drawnow;
     pause(0.1);
end

hold off;