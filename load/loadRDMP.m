%% Loads demonstrations to form rhythmic DMPs (strike & return) 
% using kinesthetic teach-in data

clc; clear; close all

% load robot and ball data and calculate striking time approx.
analyzeTableTennisDemonstrations;

for j = 1:dof
    joints{j} = ['joint\_', int2str(j)];
    vel{j} = ['joint\_vel\_', int2str(j)];
    acc{j} = ['joint\_acc\_', int2str(j)];
    %err{i} = ['err_j_', int2str(i)];
end

for i = 1:length(set)
    
    t = demo(i).t_rdmp;
    q = demo(i).Q_rdmp(1:dof,:);
    qd = demo(i).Q_rdmp(dof+1:end,:);

    % Filter the signals
    cutoff = 45;
    [q,qd,qdd] = filterDMP(t,q,qd,cutoff);

    %{
    figure;
    for j = 1:dof
        subplot(7,2,2*j-1);
        plot(t,q(j,:));
        legend(joints{j});
        subplot(7,2,2*j);
        plot(t,qd(j,:));
        legend(vel{j});
    end
    %}
    
    %close all;
    Q_rdmp{i} = q;
    Qd_rdmp{i} = qd;
    Qdd_rdmp{i} = qdd;
    t_rdmp{i} = t - t(1);


end

% feed them all to a dmp
rdmp = trainMultiDMPs(t_rdmp,Q_rdmp,Qd_rdmp,Qdd_rdmp,'r');

%% Evolve DMP and write weights to a text file

% scale up to 500 Hz 
scale = 500/200;
tf = t_rdmp{end};
tLin = linspace(tf(1),tf(end),scale*length(tf)); 
rdmp(1).can.changeSamplingTime(0.002);

Qdmp = zeros(2*dof,length(tLin));

figure('Name','RDMP Evolution');
for j = 1:dof    
    [~,qs] = rdmp(j).evolve(length(tLin));
    Qdmp(j,:) = qs(1,:);
    Qdmp(j+dof,:) = qs(2,:);
    subplot(7,3,3*j-2);
    plot(tLin,qs(1,:));
    
    legend([joints{j},'\_rdmp']);
    subplot(7,3,3*j-1);
    plot(tLin,qs(2,:));
    
    legend([vel{j},'\_rdmp']);
    subplot(7,3,3*j);
    plot(tLin,qs(3,:));
    legend([acc{j},'\_rdmp']);
end

figure('Name','RDMP in Cartesian space');
x_dmp = wam.kinematics(Qdmp);
plot3(x_dmp(1,:),x_dmp(2,:),x_dmp(3,:),'--r');
grid on;
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
legend('RDMP');

for i = 1:dof,
    W(:,i) = rdmp(i).w;
end
dlmwrite('w_rdmp.txt',W,'delimiter','\t','precision',6);

%% See how well we extend to particular demonstrations

% k = 12;
% q_act = interp1(t_rdmp{k},Q_rdmp{k}',tLin,'linear','extrap');
% qd_act = interp1(t_rdmp{k},Qd_rdmp{k}',tLin,'linear','extrap');
% qdd_act = interp1(t_rdmp{k},Qdd_rdmp{k}',tLin,'linear','extrap');
% 
% Qdmp = zeros(2*dof,length(tLin));
% Qexample = [q_act,qd_act]';
% 
% figure('Name', 'RDMP and Demonstration');
% for j = 1:dof
%     rdmp(j).setGoal(q_act(:,j));
%     rdmp(j).setInitState([q_act(1,j);0;0]);
%     rdmp(j).resetStates();
%     [~,qs] = rdmp(j).evolve(length(tLin));
%     Qdmp(j,:) = qs(1,:);
%     Qdmp(j+dof,:) = qs(2,:);
%     subplot(7,3,3*j-2);
%     plot(tLin,qs(1,:),tLin,q_act(:,j)');
%     legend([joints{j},'\_rdmp'],'demonstration');
%     subplot(7,3,3*j-1);
%     plot(tLin,qs(2,:),tLin,qd_act(:,j)');
%     legend([vel{j},'\_rdmp'],'demonstration');
%     subplot(7,3,3*j);
%     plot(tLin,qs(3,:),tLin,qdd_act(:,j)');
%     legend([acc{j},'\_rdmp'],'demonstration');
% end

%% Plot DMP and demonstrations in Cartesian space

% Qdmp = zeros(2*dof,2*length(tLin));
% 
% figure('Name','RDMP Evolved for longer');
% for j = 1:dof
%     rdmp(j).resetStates();
%     [~,qs] = rdmp(j).evolve(2*length(tLin));
%     Qdmp(j,:) = qs(1,:);
%     Qdmp(j+dof,:) = qs(2,:);
%     subplot(7,3,3*j-2);
%     plot(qs(1,:));
%     subplot(7,3,3*j-1);
%     plot(qs(2,:));
%     subplot(7,3,3*j);
%     plot(qs(3,:));
% end
% 
% x_dmp = wam.kinematics(Qdmp);
% x_example = wam.kinematics(Qexample);
% 
% figure('Name','RDMP and Demonstration in Cartesian space')
% plot3(x_dmp(1,:),x_dmp(2,:),x_dmp(3,:),'--r');
% hold on;
% plot3(x_example(1,:),x_example(2,:),x_example(3,:),'-b');
% grid on;
% axis equal;
% xlabel('x');ylabel('y');zlabel('z');
% legend('dmp',['demo ',num2str(j)]);