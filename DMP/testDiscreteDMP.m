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

%% test unforced discrete DMP

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
% type of pattern to be generated
pat = 'd';
can = Canonical(h,ax,tau,nst,pat);

% create two different DMPs
alpha = 25;
beta = 25/4;
goal = 1;
yin = [0;0];
force.w = zeros(10,1);
force.h = ones(10,1);
force.c = linspace(0,1,10);
dmp1 = discreteDMP(can,alpha,beta,goal,yin,force);
alpha = 25;
beta = 25/4;
goal = 5;
yin = [1;0];
force.w = zeros(10,1);
force.h = ones(10,1);
force.c = linspace(0,1,10);
dmp2 = discreteDMP(can,alpha,beta,goal,yin,force);

[x,y1] = dmp1.evolve();
[~,y2] = dmp2.evolve();

% plotting each dmp trajectory
t = tspan(1):h:tspan(2);
yd = y1(1,:);
y = y1(2,:);
figure(1);
plot(t,yd,'-',t,y,'-.',t,x,t,dmp1.goal*ones(1,length(t)),'r.');
legend('state velocity ydot','state y','phase','goal state');

yd = y2(1,:);
y = y2(2,:);
figure(2);
plot(t,yd,'-',t,y,'-.',t,x,t,dmp2.goal*ones(1,length(t)),'r.');
legend('state velocity ydot','state y','phase','goal state');


