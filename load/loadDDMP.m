%% Loads demonstrations to form discrete DMPs (strike only)
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
    
    t = demo(i).t_strike;
    q = demo(i).Q_strike(1:dof,:);
    qd = demo(i).Q_strike(dof+1:end,:);    

    % Compare qd recording with the numerical derivative of q
    %{
    qd2 = diff(q')./(repmat(diff(t),1,dof));
    qd2 = [qd2'; qd2(:,end)];
        
    figure;
    for j = 1:dof
        subplot(7,2,2*j-1);
        plot(t,q(:,j));
        legend(joints{j});
        subplot(7,2,2*j);
        plot(t,qd(:,j));
        legend(vel{j});
    end
    %}

    % check whether velocity data matches to differences in one joint
    %{
    figure;
    for j = 1:dof
        subplot(7,1,j);
        plot(t,qd(:,j),t,qd2(:,j));
        legend([vel{j},'\_real'],[vel{j},'\_diff']);
    end
    %}

    %% Filter the signals

    cutoff = 5;
    [qFil,qdFil,qddFil] = filterDMP(t,q,qd,cutoff);
    
    % segment motion based on maximum velocity
    % does not work due to artifacts!
%     velMax = max(abs(qdFil'));
%     cutOffVel = 0.05;
%     velThres = cutOffVel * velMax;
%     velThres = repmat(velThres(:),1,length(t));
%     idxR = abs(qd) > velThres;
%     % find the joint which lasts the longest
%     [~,idxJoint] = max(sum(idxR));
%     tNew = t(idxR(:,idxJoint));
%     qNew = qFil(idxR(:,idxJoint),:);
%     qdNew = qdFil(idxR(:,idxJoint),:);

    %{
    figure;
    for j = 1:dof
        subplot(7,3,3*j-2);
        plot(t,qFil(:,j),t,q(:,j));
        legend([joints{j},'\_filt'], joints{j})
        subplot(7,3,3*j-1);
        plot(t,qdFil(:,j),t,qd(:,j));
        legend([vel{j},'\_filt'], vel{j});
        subplot(7,3,3*j);
        plot(t,qdd(:,j));
        legend(acc{j});
    end
    %}

    %close all;
    Q_strike{i} = qFil;
    Qd_strike{i} = qdFil;
    Qdd_strike{i} = qddFil;
    t_strike{i} = t - t(1);


end

% feed them all to a dmp
dmpStrike = trainMultiDMPs(t_strike,Q_strike,Qd_strike,Qdd_strike,'d');

%% Write to text file

% evolve dmps and save to a text file
% put into matrix form
% t q1 qd1 q2 qd2 ... q7 qd7
Ms = [];

% scale up to 500 Hz 
scale = 500/200;
tf = t_strike{end};
tLin = linspace(tf(1),tf(end),scale*length(tf)); 
dmpStrike(1).can.changeSamplingTime(0.002);
Ms = tLin(:);

figure;
for j = 1:dof
    [~,qs] = dmpStrike(j).evolve(length(tLin));
    subplot(7,3,3*j-2);
    plot(tLin,qs(1,:));
    Ms = [Ms, qs(1,:)'];
    legend([joints{j},'\_dmpStrike']);
    subplot(7,3,3*j-1);
    plot(tLin,qs(2,:));
    Ms = [Ms, qs(2,:)'];
    legend([vel{j},'\_dmpStrike']);
    subplot(7,3,3*j);
    plot(tLin,qs(3,:));
    legend([acc{j},'\_dmpStrike']);
end


dlmwrite('dmp_strike.txt',Ms,'delimiter','\t','precision',6);
save('dmpStrike.mat','dmpStrike');

for i = 1:dof,
    Ws(:,i) = dmpStrike(i).w;
end
dlmwrite('w_strike.txt',Ws,'delimiter','\t','precision',6);

%% See how well we extend to particular demonstrations

k = 25;
q_act = interp1(t_strike{k},Q_strike{k}',tLin,'linear','extrap');
qd_act = interp1(t_strike{k},Qd_strike{k}',tLin,'linear','extrap');
qdd_act = interp1(t_strike{k},Qdd_strike{k}',tLin,'linear','extrap');

Qdmp = zeros(2*dof,length(tLin));
Qexample = [q_act,qd_act]';

figure;
for j = 1:dof
    dmpStrike(j).setGoal(q_act(:,j));
    dmpStrike(j).setInitState([q_act(1,j);0;0]);
    dmpStrike(j).resetStates();
    [~,qs] = dmpStrike(j).evolve(length(tLin));
    Qdmp(j,:) = qs(1,:);
    Qdmp(j+dof,:) = qs(2,:);
    subplot(7,3,3*j-2);
    plot(tLin,qs(1,:),tLin,q_act(:,j)');
    legend([joints{j},'\_dmpStrike'],'demonstration');
    subplot(7,3,3*j-1);
    plot(tLin,qs(2,:),tLin,qd_act(:,j)');
    legend([vel{j},'\_dmpStrike'],'demonstration');
    subplot(7,3,3*j);
    plot(tLin,qs(3,:),tLin,qdd_act(:,j)');
    legend([acc{j},'\_dmpStrike'],'demonstration');
end

%% Plot DMP and demonstrations in Cartesian space

x_dmp = wam.kinematics(Qdmp);
x_example = wam.kinematics(Qexample);

figure;
plot3(x_dmp(1,:),x_dmp(2,:),x_dmp(3,:),'--r');
hold on;
plot3(x_example(1,:),x_example(2,:),x_example(3,:),'-b');
grid on;
axis equal;
xlabel('x');ylabel('y');zlabel('z');
legend('dmp',['demo ',num2str(j)]);