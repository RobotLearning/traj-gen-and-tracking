%% Estimate the variance of joint positions and velocities based on kinesthetic teachin data

dof = 7; % seven degrees of freedom 
scale = 1e-3; % recorded in milliseconds

i = 15;
M = dlmread(['../Desktop/okanKinestheticTeachin_20141210/unifyData', ...
                 int2str(set(i)),'.txt']);
             
 % extract time and the joints - last 14 variables
Mq = M(:,end-14:end);

t = Mq(:,1);
t = scale * t;
q = Mq(:,2:dof+1);
qd = Mq(:,dof+2:end);

for j = 1:dof
    joints{j} = ['joint\_', int2str(j)];
    vel{j} = ['joint\_vel\_', int2str(j)];
    noise_vel{j} = ['noise_vel_j_', int2str(i)];
end

qd2 = diff(q)./(repmat(diff(t),1,dof));
qd2 = [qd2; qd2(end,:)];

figure;
for j = 1:dof
    subplot(7,2,2*j-1);
    plot(t,q(:,j));
    legend(joints{j});
    subplot(7,2,2*j);
    plot(t,qd(:,j));
    legend(vel{j});
end

% check whether velocity data matches to differences in one joint
figure;
for j = 1:dof
    subplot(7,1,j);
    plot(t,qd(:,j),t,qd2(:,j));
    legend([vel{j},'\_real'],[vel{j},'\_diff']);
end

% look at the velocity noise
for j = 1:dof
    subplot(7,1,j);
    plot(t,qd(:,j)-qd2(:,j));
    legend(noise_vel{j});
end