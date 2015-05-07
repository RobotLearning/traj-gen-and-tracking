%% Train multi-dmps

function dmp = trainMultiDMPs(t,q,qd)

% scaling up to 500 Hz
scale = 500/200;

% number of demonstrations
D = length(q);
dof = size(q{1},2);
tf = t{end};

% TODO: how to regress on nonuniform grid?
tLin = linspace(tf(1),tf(end),length(tf)); 
qLin = cell(D,1);
qdLin = cell(D,1);
qdd = cell(D,1);
goals = zeros(D,dof);

for i = 1:D

    qLin{i} = interp1(t{i},q{i},tLin,'linear','extrap');
    qdLin{i} = interp1(t{i},qd{i},tLin,'linear','extrap');

    % take numerical derivative of velocities first
    qdd{i} = diff(qdLin{i})./(repmat(diff(tLin)',1,dof));
    qdd{i}(end+1,:) = qdd{i}(end,:);
    
    goals(i,:) = qLin{i}(end,:);
end

% canonical system
h = 0.005; % 200 Hz recordings
tau = 1; %0.5/tf(end);
ax = 2; % 25/3
% number of basis functions
numbf = 50;
pat = 'd';
can = CAN(h,ax,tau,numbf,tf(end),pat);

alpha = 15;
beta = 15/4;

% stack cells 
q = cell2mat(qLin);
qd = cell2mat(qdLin);
qdd = cell2mat(qdd);

for i = 1:dof
    
    % goal and amplitude are initialized here
    goal = goals(:,i);
    % this is not important as dmps can be extended to any yin
    yin = q(1,i);
    % initial states of DMPs
    dmp(i) = DDMP(can,alpha,beta,goal(1),yin);
    dmp(i).regressLive(q(:,i),qd(:,i),qdd(:,i),goal);
    
end

% scale evolution up to 500 Hz
dmp(1).can.dt = h/scale; 