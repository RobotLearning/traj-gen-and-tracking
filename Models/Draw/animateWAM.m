%% creates an animation for the Barrett WAM

function animateWAM(x,Rracket)

% x1 and x2 are the cartesian coordinates of joint 2 and end-effector
link1_x = [shift(1) x1(1,1)];
link1_y = [shift(2) x1(2,1)];
link2_x = [x1(1,1),x2(1,1)];
link2_y = [x1(2,1),x2(2,1)];
hf = figure('color','white');
axis manual; axis equal; 

h1 = line(link1_x, link1_y, 'color', [.4 .4 .4],'LineWidth',4);
h2 = line(link2_x, link2_y, 'color', [.4 .4 .4],'LineWidth',4);
title('Two-link planar arm');
hold on;
h3 = scatter(shift(1),shift(2),100,'b','LineWidth',4);
h4 = scatter(x1(1,1),x1(2,1),100,'b','LineWidth',4);
h5 = scatter(x2(1,1),x2(2,1),10,'b','LineWidth',4);
h6 = plot(x_des,y_des,'-.','color', [.6 .6 .6], 'LineWidth',1);

for j = 2:length(x)
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
     %traj_x = [x2(1,j-1),x2(1,j)];               
     %traj_y = [x2(2,j-1),x2(2,j)];              
     %h7 = line(traj_x, traj_y, 'color', [.4 .4 .5],'LineWidth',1);
     
     drawnow;
     pause(0.005);
end

%hold off;