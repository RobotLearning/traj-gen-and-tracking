%% Train multi-dmps

function dmp = trainMultiDMPs(t,q,qd,qdd)

f = 200; %200 Hz recording

% number of demonstrations
D = length(q);
dof = size(q{1},2);
tf = t{end};

% TODO: how to regress on nonuniform grid?
tLin = linspace(tf(1),tf(end),length(tf)); 
qLin = cell(D,1);
qdLin = cell(D,1);
qddLin = cell(D,1);
goals = zeros(D,dof);
inits = zeros(D,dof);

for i = 1:D

    qLin{i} = interp1(t{i},q{i},tLin,'linear','extrap');
    qdLin{i} = interp1(t{i},qd{i},tLin,'linear','extrap');
    qddLin{i} = interp1(t{i},qdd{i},tLin,'linear','extrap');

    % numerical differentiation of velocities first
    %qdd{i} = diff(qdLin{i})./(repmat(diff(tLin)',1,dof));
    %qdd{i}(end+1,:) = qdd{i}(end,:);
    
    goals(i,:) = qLin{i}(end,:);
    vels(i,:) = qdLin{i}(end,:);
    accs(i,:) = qddLin{i}(end,:);
    inits(i,:) = qLin{i}(1,:);
end

% canonical system
h = 1/f; % 200 Hz recordings
tau = 1/tf(end);
alpha = 15;
beta = alpha/4;
ax = 1;
% number of basis functions
numbf = 50;
pat = 'd';
can = CAN(h,ax,tau,numbf,tf(end),pat);

% stack cells 
q = cell2mat(qLin);
qd = cell2mat(qdLin);
qdd = cell2mat(qddLin);

g = zeros(1,dof);
yin = zeros(3,dof);

for i = 1:dof
    
    % tau is fixed to be one
    % goal and amplitude are initialized here
    % these are not important as dmps can be extended to any yin and goal
    yin(:,i) = [inits(1,i);0;0];
    g(i) = goals(1,i);
    % initial states of DMPs
    dmp(i) = DDMP(can,alpha,beta,g(i),yin(:,i));
    dmp(i).regressLive(q(:,i),qd(:,i),qdd(:,i),goals(:,i));
    
end

disp('Initializing DMP to yin = ');
disp(yin(1,:)');
disp('Goal state for DMP, goal = ');
disp(g(:));