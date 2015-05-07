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
tspan = [0 1];
% number of points
t = tspan(1):h:tspan(2);
N = length(t);
% number of basis functions
numbf = 50;
% type of pattern to be generated
pat = 'r';
can = CAN(h,ax,tau,numbf,tspan(2)-tspan(1),pat);

%% test forced rhythmic DMP

% create two different DMPs
alpha = 25;
beta = 25/4;
% goal and amplitude are initialized here
goal = 1;
amp = 1;
% initial states of DMPs
yin1 = 0;
yin2 = 0;
dmp1 = RDMP(can,alpha,beta,goal,amp,yin1);
dmp2 = RDMP(can,alpha,beta,goal,amp,yin2);

% create two paths
ctr1 = 5; % center 
amp1 = 2; % amplitude
path1 = ctr1 + amp1 * sin(10*pi*t);
ctr2 = 1;
amp2 = 2;
path2 = ctr2 + amp2 * cos(10*pi*t);

% learn the weights with the usual linear regression
dmp1.LWR(path1(end),path1);
dmp2.LWR(path1(end),path2);

[x,y1] = dmp1.evolve(N);
[~,y2] = dmp2.evolve(N);

% plotting each dmp trajectory
y = y1(1,:);
figure(3);
plot(t,path1,'-',t,y,'-.',t,sin(x));
legend('desired trajectory','state y','sin of phase');
title('Followed trajectory for DMP1');

y = y2(1,:);
figure(4);
plot(t,path2,'-',t,y,'-.',t,sin(x));
legend('desired trajectory','state y','sin of phase');
title('Followed trajectory for DMP2');