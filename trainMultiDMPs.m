%% Train multi-dmps

function dmp = trainMultiDMPs(t,q,qd)

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
h = tLin(2)-tLin(1);
tau = 1;
ax = 2;
nst = length(tf);
pat = 'd';
can = Canonical(h,ax,tau,nst,pat);

alpha = 25;
beta = 25/4;
% number of basis functions
numbf = 50;

% stack cells 
q = cell2mat(qLin);
qd = cell2mat(qdLin);
qdd = cell2mat(qdd);

for i = 1:dof
    
    % goal and amplitude are initialized here
    goal = goals(:,i);
    % this is not important as dmps can be extended to any yin
    yin = [q(1,i);qd(1,i)];
    % initial states of DMPs
    dmp(i) = discreteDMP(can,alpha,beta,goal(1),yin,numbf);
    dmp(i).regressLive(q(:,i),qd(:,i),qdd(:,i),goal);
    
end