%% creates an animation for the two-wheeled robot kinematical model

function animateCar(X,s,PAR)

% extract parameters
scale = 0.1;
d = scale * PAR.length;
r1 = scale * PAR.wheel1.radius;
r2 = scale * PAR.wheel2.radius;

% extract desired trajectory
x_des = s(1,:);
y_des = s(2,:);

% extract followed trajectory
x = X(1,:);
y = X(2,:);
theta = X(3,:);

% make animations of two wheels and one connecting link

link_x = [x(1) - d*sin(theta(1)), x(1) + d*sin(theta(1))];
link_y = [y(1) + d*cos(theta(1)), y(1) - d*cos(theta(1))];
wheel1_x = [link_x(1) - r1*cos(theta(1)), link_x(1) + r1*cos(theta(1))];
wheel1_y = [link_y(1) - r1*sin(theta(1)), link_y(1) + r1*sin(theta(1))];
wheel2_x = [link_x(2) - r2*cos(theta(1)), link_x(2) + r2*cos(theta(1))];
wheel2_y = [link_y(2) - r2*sin(theta(1)), link_y(2) + r2*sin(theta(1))];

hf = figure('color','white');
%axis manual, axis equal;
h1 = line(link_x, link_y, 'color', [.4 .4 .8],'LineWidth',4);
h2 = line(wheel1_x, wheel1_y, 'color', [.4 .4 .8],'LineWidth',4);
h3 = line(wheel2_x, wheel2_y, 'color', [.4 .4 .8],'LineWidth',4);
title('Two-wheeled robot kinematics model','Color',[.6 0 0]);
hold on;
h4 = plot(x_des,y_des,'-r','LineWidth',2);

for j = 2:length(s)
     link_x(1) = x(j) - d*sin(theta(j));
     link_x(2) = x(j) + d*sin(theta(j));
     link_y(1) = y(j) + d*cos(theta(j));
     link_y(2) = y(j) - d*cos(theta(j));
     wheel1_x(1) = link_x(1) - r1*cos(theta(j));
     wheel1_x(2) = link_x(1) + r1*cos(theta(j));
     wheel1_y(1) = link_y(1) - r1*sin(theta(j));
     wheel1_y(2) = link_y(1) + r1*sin(theta(j));
     wheel2_x(1) = link_x(2) - r2*cos(theta(j));
     wheel2_x(2) = link_x(2) + r2*cos(theta(j));
     wheel2_y(1) = link_y(2) - r2*sin(theta(j));
     wheel2_y(2) = link_y(2) + r2*sin(theta(j));
     
     set(h1,'XData',link_x)
     set(h1,'YData',link_y)
     set(h2,'XData',wheel1_x)
     set(h2,'YData',wheel1_y)
     set(h3,'XData',wheel2_x);
     set(h3,'YData',wheel2_y);
     
     % connect the resulting end-effector trajectory
     traj_x = [x(j-1),x(j)];               
     traj_y = [y(j-1),y(j)];              
     h5 = line(traj_x, traj_y, 'color', [.4 .4 .5],'LineWidth',2);
     
     drawnow;
     pause(0.01);
end

hold off;