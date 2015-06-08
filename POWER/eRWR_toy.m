% This file contains a sample implementation of a episodic
% Reward Weighted Regression based policy search method
% The required motor primitive code can be downloaded from
% http://www-clmc.usc.edu/Resources/Software

clear all;
close all;
clc;

% time step
dt = 0.005;

% length of the episode
T = 1;

% number of basis functions
n_rfs = 10;

ID = 1;
goal = 1;

% number of iterations
n_iter = 1000;

% initialize the motor primitive and storage variables
dcp('init',ID,n_rfs,'eRWR');

rt = 0:dt:T;

y = zeros(1,length(rt));
yd = zeros(1,length(rt));
ydd = zeros(1,length(rt));

explore = zeros(1,length(rt));
basis = zeros(n_rfs,length(rt));

Return = zeros(1,n_iter+1);
s_Return = zeros(n_iter+1,2);
param = zeros(n_rfs,n_iter+1);

pre_dnom=zeros(n_rfs,n_iter);
pre_nom=zeros(n_rfs,n_iter);

randn('state',20);

% set the standard deviation for the exploration
std_dev = 50.;

% start with all parameters set to zero
param(:,1) = zeros(n_rfs,1); 

% apply the new parameters to the motor primitive
dcp('change',ID,'w',param(:,1));

% reset the motor primitive
dcp('reset_state',ID);
dcp('set_goal',ID,goal,1);

% run the motor primitive & precalculate the basis functions
for i=1:length(rt)
    [y(i),yd(i),ydd(i),basis(:,i)] = dcp('run',ID,T,dt);
end
% dcp('run') returns the unnormalized basis functions
basis = basis./(ones(n_rfs,1)*sum(basis,1));

% precalculate the canonical system
x = zeros(1,length(rt));
x(1) = 1;
for i=2:length(rt),
    xd   = -25/3*x(i-1)*(.5/T);
    x(i) = x(i-1)+xd*dt;
end

% do the iterations
for iter=1:n_iter
    if (mod(iter,100)==0)
        disp(['Iter ', num2str(iter)]);
    end
    
    % calculate the return of the previous rollout
    % the policy starts at 0 and the distance to 1 is punished
    Return(iter) = sum(exp(-abs(1-y)))/length(rt);
    
    % precompute for the reward weighted regression
    fT = (ydd/(.5/T)^2-25*(25/4*(goal-y)-yd/(.5/T)));
    pre_dnom(:,iter) = Return(iter)*sum(basis.*(ones(n_rfs,1)*(x.^2)),2);
    pre_nom(:,iter) = Return(iter)*sum(basis.*(ones(n_rfs,1)*(x.*fT)),2);
    
    % this lookup table will be used for the importance sampling
    s_Return(1,:) = [Return(iter) iter];
    s_Return = sortrows(s_Return);
    
    % update the policy parameters
    param_nom = zeros(n_rfs,1);
    param_dnom = zeros(n_rfs,1);
    
    % as importance sampling we take the 10 best rollouts
    for i=1:min(iter,10)
        % get the rollout number for the 10 best rollouts
        j = s_Return(end+1-i,2);
        
        param_nom=param_nom+pre_nom(:,j);
        param_dnom=param_dnom+pre_dnom(:,j);
    end
    
    % update the parameters
    param(:,iter+1) = param_nom./(param_dnom+1.e-10);
    
    % generate the exploration
    explore(:) = (randn(1,length(rt))).*std_dev;
    
    % in the last rollout we want to get the return without exploration
    if iter==n_iter
        explore(:) = zeros(1,length(rt));
    end
    
    % apply the new parameters to the motor primitve
    dcp('change',ID,'w',param(:,iter+1));
    
    % reset the motor primitive
    dcp('reset_state',ID);
    dcp('set_goal',ID,goal,1);
    
    % run the motor primitive
    % run the motor primitive
    for i=1:length(rt)
        [y(i),yd(i),ydd(i)] = dcp('run',ID,T,dt,explore(i));
    end
    
end

% calculate the return of the final rollout
Return(iter+1) = sum(exp(-abs(1-y)))/length(rt);

% plot the return over the rollouts
figure(1);
plot(Return);
ylabel('return');
xlabel('rollouts');

disp(['Final Return ', num2str(Return(end))]);
