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
ax = 2;
% temporal scaling factor
tau = 1;
% time evolution
tspan = [0 1];
% number of points
nst = tspan(2)/h + 1;
t = tspan(1):h:tspan(2);
% type of pattern to be generated
pat = 'd';
can = Canonical(h,ax,tau,nst,pat);

%% test unforced discrete DMP

% create two different DMPs
alpha = 25;
beta = 25/4;
goal = 1;
yin = [0;0];
force.w = zeros(10,1);
force.h = ones(10,1);
force.c = linspace(0,1,10);
dmp1 = discreteDMP(can,alpha,beta,goal,yin,force);

goal = 5;
yin = [1;0];
force.w = zeros(10,1);
force.h = ones(10,1);
force.c = linspace(0,1,10);
dmp2 = discreteDMP(can,alpha,beta,goal,yin,force);

[x,y1] = dmp1.evolve();
[~,y2] = dmp2.evolve();

% plotting each dmp trajectory
y = y1(1,:);
yd = y1(2,:);
figure(1);
plot(t,yd,'-',t,y,'-.',t,x,t,dmp1.goal*ones(1,length(t)),'r.');
legend('state velocity ydot','state y','phase','goal state');
title('Unforced trajectory for DMP1');

y = y2(1,:);
yd = y2(2,:);
figure(2);
plot(t,yd,'-',t,y,'-.',t,x,t,dmp2.goal*ones(1,length(t)),'r.');
legend('state velocity ydot','state y','phase','goal state');
title('Unforced trajectory for DMP2');

%% test forced discrete DMP

% create two paths
goal1 = 5;
path1 = goal1*t.*sin(5*pi/2 * t);
goal2 = 1;
jump2 = floor(length(t)/2);
path2 = [zeros(1,jump2), goal2*ones(1,length(t)-jump2)];

alpha = 25;
beta = 25/4;

% learn the weights with locally weighted regression
force.h = ones(100,1) * 100^(1.5);
force.c = linspace(0,1,100);
force1 = LWR(path1,can,alpha,beta,force);
force2 = LWR(path2,can,alpha,beta,force);

% create two different DMPs
yin1 = [path1(1);0];
dmp1 = discreteDMP(can,alpha,beta,goal1,yin1,force1);
yin2 = [path2(1);0];
dmp2 = discreteDMP(can,alpha,beta,goal2,yin2,force2);

[x,y1] = dmp1.evolve();
[~,y2] = dmp2.evolve();

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