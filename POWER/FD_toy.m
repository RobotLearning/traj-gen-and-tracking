% This file contains a sample implementation of a Finite Difference
% Gradient based policy search method
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
% a lot more than for PoWER...
n_iter = 5000;

% batch size
n_batch = n_rfs+1;
n_iter = floor(n_iter/n_batch);

% step size
step = 100000;

% initialize the motor primitive and storage variables
dcp('init',ID,n_rfs,'FD');

rt = 0:dt:T;

y = zeros(1,length(rt));
yd = zeros(1,length(rt));
ydd = zeros(1,length(rt));

explore = zeros(n_rfs,n_batch);
batch_return = zeros(n_batch,1);
Return = zeros(1,n_iter*n_batch+2);

param = zeros(n_rfs,n_iter+1);
rand('twister',20);

% set the exploration magnitude
explore_magn = .01*ones(n_rfs,1);

% start with all parameters set to zero
param(:,1) = zeros(n_rfs,1); 

% apply the new parameters to the motor primitive
dcp('change',ID,'w',param(:,1));

% reset the motor primitive
dcp('reset_state',ID);
dcp('set_goal',ID,goal,1);

% run the motor primitive
for i=1:length(rt)
    [y(i),yd(i),ydd(i)] = dcp('run',ID,T,dt);
end

% calculate the return of the rollout
Return(1) = sum(exp(-abs(1-y)))/length(rt);

% do the iterations
for iter=1:n_iter
    if (mod(iter,10)==1 && iter>1)
        disp(['Iter ', num2str((iter-1)*n_batch)]);
    end
    
    for batch=1:n_batch
        % generate the exploration
        explore(:,batch) = (rand(n_rfs,1)-.5).*explore_magn;
        % apply the new parameters to the motor primitve
        dcp('change',ID,'w',param(:,iter)+explore(:,batch));
        % reset the motor primitive
        dcp('reset_state',ID);
        dcp('set_goal',ID,goal,1);
        
        % run the motor primitive
        for i=1:length(rt)
            [y(i),yd(i),ydd(i)] = dcp('run',ID,T,dt);
        end
        
        % calculate the return of the batch rollout
        % the policy starts at 0 and the distance to 1 is punished
        batch_return(batch) = sum(exp(-abs(1-y)))/length(rt);
        Return((iter-1)*n_batch+batch+1) = sum(exp(-abs(1-y)))/length(rt);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    explore_temp=[explore; ones(1,batch)];
    gradient=inv(explore_temp*explore_temp'+eye(batch)*1.e-10)*explore_temp*batch_return;
    param(:,iter+1) = param(:,iter)+step*gradient(1:10);
    
    % apply the new parameters to the motor primitve
    dcp('change',ID,'w',param(:,iter+1));    
end

% reset the motor primitive
dcp('reset_state',ID);
dcp('set_goal',ID,goal,1);

% run the motor primitive
for i=1:length(rt)
    [y(i),yd(i),ydd(i)] = dcp('run',ID,T,dt);
end

% calculate the return of the rollout
Return(end) = sum(exp(-abs(1-y)))/length(rt);

% plot the return over the rollouts
figure(1);
plot(Return);
ylabel('return');
xlabel('rollouts');

disp(['Final Return ', num2str(Return(end))]);
