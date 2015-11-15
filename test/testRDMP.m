%% test script for rhythmic DMP

%# store breakpoints
tmp = dbstatus;
save('tmp.mat','tmp')

%# clear all
close all
clear classes %# clears even more than clear all
clc

%# reload breakpoints
load('tmp.mat')
dbstop(tmp)

%# clean up
clear tmp
delete('tmp.mat')

%% create canonical system

% discretization value - dt
h = 0.01; 
% time constant of the canonical system
% this decay rate has no effect on the rhythmic system
ax = 1;
% temporal scaling factor
tau = 1;
% time evolution
T = 3;
tspan = [0 T];
% number of points
t = tspan(1):h:tspan(2);
N = length(t);
% number of basis functions
numbf = 50;
% type of pattern to be generated
pat = 'r';
can = CAN(h,ax,tau,numbf,T,pat);

%% test forced rhythmic DMP

% create two different DMPs
alpha = 50;
beta = 50/4;
% goal and amplitude are initialized here
goal = 1;
amp = 1;
% initial states of DMPs
yin1 = [5;0;0];
yin2 = [16;0;0];
dmp1 = RDMP(can,alpha,beta,goal,amp,yin1);
dmp2 = RDMP(can,alpha,beta,goal,amp,yin2);

% create two paths
ctr1 = 5; % center 
amp1 = 5; % amplitude
path1 = ctr1 + amp1 * sin(2*pi*t);
ctr2 = 7;
amp2 = 10;
path2 = ctr2 + amp2 * cos(2*pi*t);

% learn the weights with the usual linear regression
dmp1.setGoal(path1);
dmp2.setGoal(path2);
dmp1.updateWeights(path1);
dmp2.updateWeights(path2);

[x,y1] = dmp1.evolve(N);
[~,y2] = dmp2.evolve(N);

% plotting each dmp trajectory
y1 = y1(1,:);
figure(1);
plot(t,path1,'-',t,y1,'-.',t,sin(x));
legend('desired trajectory','state y','sin of phase');
title('Followed trajectory for DMP1');

y2 = y2(1,:);
figure(2);
plot(t,path2,'-',t,y2,'-.',t,sin(x));
legend('desired trajectory','state y','sin of phase');
title('Followed trajectory for DMP2');

figure(3);
plot(path1,path2,'-',y1,y2,'-.');
legend('path','dmp');

%% test regression on Barrett WAM

% initializeWAM;
% q0 = zeros(N_DOFS,3);
% goal = ones(N_DOFS,1);
% amp = ones(N_DOFS,1);
% 
% for j = 1:N_DOFS
%     joints{j} = ['joint\_', int2str(j)];
%     vel{j} = ['joint\_vel\_', int2str(j)];
%     acc{j} = ['joint\_acc\_', int2str(j)];
%     %err{i} = ['err_j_', int2str(i)];
% end
% 
% q = zeros(N_DOFS,N);
% qd = zeros(N_DOFS,N);
% qdd = zeros(N_DOFS,N);
% 
% for i = 1:N_DOFS
%     dmp(i) = RDMP(can,alpha,beta,goal,amp,q0(i,:)');
%     q(i,:) = rand + rand * cos(2*pi*t) + rand * sin(2*pi*t);
%     qd(i,1:end-1) = diff(q(i,:))/h; qd(i,end) = qd(i,end-1);
%     qdd(i,1:end-1) = diff(qd(i,:))/h; qdd(i,end) = qdd(i,end-1);
%     dmp(i).regressLive(q(i,:)',qd(i,:)',qdd(i,:)',goal(i));
% end
% 
% figure;
% for i = 1:N_DOFS
%     dmp(i).resetStates();
%     [~,qdmp] = dmp(i).evolve(N);
%     Qdmp(i,:) = qdmp(1,:);
%     Qdmp(i+N_DOFS,:) = qdmp(2,:);
%     subplot(7,2,2*i-1);
%     plot(t,qdmp(1,:),t,q(i,:));
%     legend([joints{i},'\_rdmp'],'demonstration');
%     subplot(7,2,2*i);
%     plot(t,qdmp(2,:),t,qd(i,:));
%     legend([vel{i},'\_rdmp'],'demonstration');
% 
% end
% 
% x_dmp = wam.kinematics(Qdmp);
% x_example = wam.kinematics([q;qd]);
% 
% figure;
% plot3(x_dmp(1,:),x_dmp(2,:),x_dmp(3,:),'--r');
% hold on;
% plot3(x_example(1,:),x_example(2,:),x_example(3,:),'-b');
% grid on;
% axis equal;
% xlabel('x');ylabel('y');zlabel('z');
% legend('dmp',['demo ',num2str(j)]);