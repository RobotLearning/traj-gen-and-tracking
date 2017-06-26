%% creates an animation for the pendulum

function animateR(x,s)


% rotate everything by 90 degrees
R = [0 1; -1 0];

% extract desired trajectory
x_des = s(1,:);
y_des = s(2,:);

link_x = [0 x(1,1)];
link_y = [0 x(2,1)];
hf = figure('color','white');
axis manual; axis equal; 

h1 = line(link_x, link_y, 'color', [.6 .6 .6],'LineWidth',4);
title('Pendulum animation');
hold on;
h2 = scatter(0,0,100,'b','LineWidth',4);
h3 = scatter(x(1,1),x(2,1),100,'b','LineWidth',4);
h4 = plot(x_des,y_des,'-.','color', [.6 .6 .6], 'LineWidth',1);

for j = 2:length(x)
    
     link_x(2) = x(1,j);
     link_y(2) = x(2,j);
     
     set(h1,'XData',link_x)
     set(h1,'YData',link_y)
     set(h3,'XData',link_x(2));
     set(h3,'YData',link_y(2));
     
     % connect the resulting end-effector trajectory
     %traj_x = [x2(1,j-1),x2(1,j)];               
     %traj_y = [x2(2,j-1),x2(2,j)];              
     %h = line(traj_x, traj_y, 'color', [.4 .4 .5],'LineWidth',1);
     
     drawnow;
     pause(0.005);
end

hold off;