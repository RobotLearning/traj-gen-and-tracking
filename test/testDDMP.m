%% test script for discrete DMP

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
ax = 25/3;
% temporal scaling factor
tau = 1;
% time evolution
tin = 0;
tf = 1;
% tspan = [tin tf];
% number of basis functions
nst = 50;
t = tin:h:tf;
N = length(t);
% type of pattern to be generated
pat = 'd';
can = CAN(h,ax,tau,nst,tf,pat);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

%% test unforced discrete DMP

% % create two different DMPs
% alpha = 25;
% beta = 25/4;
% goal = 1;
% yin = 0;
% dmp1 = DDMP(can,alpha,beta,goal,yin);
% 
% goal = 5;
% yin = 1;
% dmp2 = DDMP(can,alpha,beta,goal,yin);
% 
% [x,y1] = dmp1.evolve(N);
% [~,y2] = dmp2.evolve(N);
% 
% % plotting each dmp trajectory
% y = y1(1,:);
% figure(1);
% plot(t,y,'-.',t,x,t,dmp1.goal*ones(1,length(t)),'r.');
% legend('state y','phase','goal state');
% title('Unforced trajectory for DMP1');
% 
% y = y2(1,:);
% figure(2);
% plot(t,y,'-.',t,x,t,dmp2.goal*ones(1,length(t)),'r.');
% legend('state y','phase','goal state');
% title('Unforced trajectory for DMP2');

%% test forced discrete DMP and compare LWR vs. simple Regression

% create two different DMPs
alpha = 25;
beta = 25/4;
% goal and amplitude are initialized here
goal = 1;
% initial states of DMPs
yin = [0;0;0];
dmp1 = DDMP(can,alpha,beta,goal,yin);
dmp2 = DDMP(can,alpha,beta,goal,yin);

% create two paths
path1 = 5*t.*sin(5*pi/2 * t);
jump2 = floor(length(t)/2);
path2 = [zeros(1,jump2), ones(1,length(t)-jump2)];

% learn the weights with the usual linear regression
dmp1.LWR(goal,path1);
dmp2.LWR(goal,path2);

[x,y1] = dmp1.evolve(N);
[~,y2] = dmp2.evolve(N);

% plotting each dmp trajectory
y = y1(1,:);
figure(3);
plot(t,path1,'-',t,y,'-.',t,x);
legend('desired trajectory','state y','phase');
title('Followed trajectory for DMP1');

y = y2(1,:);
figure(4);
plot(t,path2,'-',t,y,'-.',t,x);
legend('desired trajectory','state y','phase');
title('Followed trajectory for DMP2');

%% test circular path

% create two paths
path1 = cos(2*pi*t);
path2 = sin(2*pi*t);

% learn the weights with the usual linear regression
dmp1.regression(1,path1);
dmp2.regression(0,path2);

% initial states of DMPs
%yin1 = [1;0];
%yin2 = [0;2*pi];
yin1 = [1;0;0];
yin2 = [0;0;0];
dmp1.setInitState(yin1);
dmp1.setGoal(path1);
dmp1.resetStates();
dmp2.setInitState(yin2);
dmp2.setGoal(path2);
dmp2.resetStates();

[x,y1] = dmp1.evolve(N);
[~,y2] = dmp2.evolve(N);

% plotting each dmp trajectory
y11 = y1(1,:);
figure(5);
plot(t,path1,'-',t,y11,'-.',t,x);
legend('desired trajectory','state y','phase');
title('Followed trajectory for DMP1');

y12 = y2(1,:);
figure(6);
plot(t,path2,'-',t,y12,'-.',t,x);
legend('desired trajectory','state y','phase');
title('Followed trajectory for DMP2');

% plot the x-y trajectory tracked
figure(7);
plot(path1,path2,'-',y11,y12,'-.');
legend('desired trajectory','followed trajectory');
title('Trajectory in xy space');

%% Test the linear model s = Fs * w

% close all;
% % create two paths
% path1 = cos(2*pi*t);
% path2 = sin(2*pi*t);
% 
% % goal position
% goal = path1(end);
% 
% % learn the weights with the usual linear regression
% dmp1.setWeights(path1);
% dmp2.setWeights(path2);
% 
% yin1 = [0;0];
% yin2 = [0;0];
% dmp1.setInitState(yin1);
% dmp1.resetStates();
% dmp2.setInitState(yin2);
% dmp2.resetStates();
% [~,sref] = dmp1.evolve();
% 
% % construct As
% tau = dmp1.can.tau;
% As = tau * [0, 1, 0; 
%      -alpha*beta, -alpha, alpha*beta; 
%             0, 0, 0];
%         
% % construct Phi
% dmp1.can.reset();
% N = length(t)-1;
% M = length(dmp1.FORCE.w);
% Phi = zeros(N,M);
% for i = 1:N
%     sumphi = 0;
%     for j = 1:M
%         x = dmp1.can.x;
%         Phi(i,j) = x * dmp1.basis(x,dmp1.FORCE.h(j),dmp1.FORCE.c(j));
%         sumphi = sumphi + dmp1.basis(x,dmp1.FORCE.h(j),dmp1.FORCE.c(j));
%     end
%     Phi(i,:) = Phi(i,:)/sumphi;
%     dmp1.can.step(1);
% end
% 
% % discretize As and Phi
% As = eye(3) + h * As;
% Phi = h * Phi;
% 
% s1 = [yin1;goal];
% dim_x = 2;
% dim_w = M;
% Fs = zeros(N*dim_x,N*M);
% % add the evolution of s0
% sfree = zeros(N*2,1);
% sf = s1;
% % construct Fs
% for i = 1:N
%     vec_x = (i-1)*dim_x + 1:i*dim_x;
%     for j = 1:i        
%         vec_w = (j-1)*dim_w + 1:j*dim_w;
%         % put zeros in between
%         mat = [zeros(1,M); Phi(j,:)];
%         for k = j+1:i
%             mat = As(1:2,1:2) * mat;
%         end
%         Fs(vec_x,vec_w) = mat; 
%     end
%     sf = As * sf;
%     sfree(vec_x) = sf(1:dim_x);
% end
% 
% % get the weights
% Fs = Fs * repmat(eye(M),N,1);
% s = sref(:,2:end);
% w_est = pinv(Fs) * (s(:) - sfree);
% w = dmp1.FORCE.w;
% % max diff
% max(abs(w_est - w))

%% Test changing the goal position

% create two paths
path1 = log(t+1);
path2 = sin(2*pi*t);

% learn the weights with the usual linear regression
dmp1.regression(path1(end),path1);
dmp2.regression(path2(end),path2);

% initial states of DMPs
yin1 = [0;0;0];
yin2 = [1;0;0];
dmp1.setInitState(yin1);
dmp1.setGoal(path1);
dmp1.resetStates();
dmp2.setInitState(yin2);
dmp2.setGoal(path2);
dmp2.resetStates();

[x,y1] = dmp1.evolve(N);
[~,y2] = dmp2.evolve(N);
y11 = y1(1,:);
y12 = y2(1,:);

dmp1.resetStates();
dmp2.resetStates();
dmp1.setGoal(path1(end));
dmp2.setGoal(1);

[~,y3] = dmp1.evolve(N);
[~,y4] = dmp2.evolve(N);
y21 = y3(1,:);
y22 = y4(1,:);

% plot the x-y trajectory tracked
figure(8);
plot(path1,path2,'-',y11,y12,'-',y21,y22,'-.');
legend('desired trajectory','original DMP','stretched DMP');
title('Trajectory in xy space');