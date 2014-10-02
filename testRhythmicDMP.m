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
tspan = [0 10];
% number of points
nst = tspan(2)/h + 1;
t = tspan(1):h:tspan(2);
% type of pattern to be generated
pat = 'r';
can = Canonical(h,ax,tau,nst,pat);

%% test forced rhythmic DMP

% create two different DMPs
alpha = 25;
beta = 25/4;
% number of basis functions
numbf = 50;
force.h = ones(numbf,1) * numbf^(1.5);
force.c = linspace(0,2*pi,numbf);
% goal and amplitude are initialized here
goal = 1;
amp = 1;
% initial states of DMPs
yin1 = [0;0];
yin2 = [0;0];
dmp1 = rhythmicDMP(can,alpha,beta,goal,amp,yin1,force);
dmp2 = rhythmicDMP(can,alpha,beta,goal,amp,yin2,force);

% create two paths
ctr1 = 5; % center 
amp1 = 2; % amplitude
path1 = ctr1 + amp1 * sin(pi*t);
ctr2 = 1;
amp2 = 2;
path2 = ctr2 + amp2 * cos(pi*t);

% learn the weights with locally weighted regression
dmp1 = LWR(path1,dmp1,force);
dmp2 = LWR(path2,dmp2,force);

[x,y1] = dmp1.evolve();
[~,y2] = dmp2.evolve();

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