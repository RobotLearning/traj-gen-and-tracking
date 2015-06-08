% This file contains a sample implementation of PoWER
% The required motor primitive code can be downloaded from
% http://www-clmc.usc.edu/Resources/Software
% 
% The algorithm is implemented with a number of simplifications:
% - the variance of the exploration is constant over trials
% - the exploration is constant during the trial
%   (as the motor primitives employ basis functions that are localized in time 
%   which are only active for a short period of time,
%   time-varying exploration does not have large effetcs)
% - always only one basis function is active

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
dcp('init',ID,n_rfs,'PoWER');

rt = 0:dt:T;

y = zeros(1,length(rt));
yd = zeros(1,length(rt));
ydd = zeros(1,length(rt));

Return = zeros(1,n_iter+1);
Q = zeros(length(rt),n_iter+1);
s_Return = zeros(n_iter+1,2);
param = zeros(n_rfs,n_iter+1);
basis = zeros(n_rfs,length(rt));
randn('state',20);

% set the variance
variance = 4000.*ones(n_rfs,1);

% start with all parameters set to zero
param(:,1) = zeros(n_rfs,1); 
current_param = param(:,1);

% apply the new parameters to the motor primitive
dcp('change',ID,'w',param(:,1));

% reset the motor primitive
dcp('reset_state',ID);
dcp('set_goal',ID,goal,1);

% run the motor primitive & precalculate the basis functions
for i=1:length(rt)
	% also store the values of the basis functions
    [y(i),yd(i),ydd(i),basis(:,i)] = dcp('run',ID,T,dt);
	% calculate the Q values
    Q(1:(end+1-i),1) = Q(1:(end+1-i),1) + exp(-abs(1-y(i)));
end
% normalize the Q values
Q(:,1) = Q(:,1)./length(rt);

% dcp('run') returns the unnormalized basis functions
basis = basis'./(sum(basis,1)'*ones(1,n_rfs));

% do the iterations
for iter=1:n_iter
    if (mod(iter,100)==0)
        disp(['Iter ', num2str(iter)]);
    end
    
    Return(iter) = Q(1,iter);

	% this lookup table will be used for the importance sampling
    s_Return(1,:) = [Return(iter) iter];
    s_Return = sortrows(s_Return);
    
    % update the policy parameters
    param_nom = zeros(n_rfs,1);
    param_dnom = zeros(n_rfs,1);

	% calculate the expectations (the normalization is taken care of by the division)
    % as importance sampling we take the 10 best rollouts
    for i=1:min(iter,10)
        % get the rollout number for the 10 best rollouts
        j = s_Return(end+1-i,2);
		% calulate W with the assumption that always only one basis functions is active
        temp_W = (basis.^2./(sum(basis.^2.*(ones(length(rt),1)*variance'),2)*ones(1,n_rfs)))';
		% calculate the exploration with respect to the current parameters
        % if you have time-varying exploration use 
        % temp_explore = (reshape(param(:,:,j),length(rt),n_rfs)-ones(length(rt),1)*current_param')';
        % instead
        temp_explore = (ones(length(rt),1)*(param(:,j)-current_param)')';
		% repeat the Q values
        temp_Q = (Q(:,j)*ones(1,n_rfs))';
        param_nom = param_nom + sum(temp_W.*temp_explore.*temp_Q,2);
        param_dnom = param_dnom + sum(temp_W.*temp_Q,2);
    end
    
    % update the parameters
    param(:,iter+1) = current_param + param_nom./(param_dnom+1.e-10);
    % set the new mean of the parameters
    current_param = param(:,iter+1);
    
    % in the last rollout we want to get the return without exploration
    if iter~=n_iter
        param(:,iter+1) = param(:,iter+1) + variance.^.5.*randn(n_rfs,1);
    end
    
    % apply the new parameters to the motor primitve
    dcp('change',ID,'w',param(:,iter+1));
    
    % reset the motor primitive
    dcp('reset_state',ID);
    dcp('set_goal',ID,goal,1);
    
    % run the motor primitive
    for i=1:length(rt)
        [y(i),yd(i),ydd(i)] = dcp('run',ID,T,dt);
		% calculate the Q values
        Q(1:(end+1-i),iter+1) = Q(1:(end+1-i),iter+1) + exp(-(1-y(i)).^2);
    end
	% normalize the Q values
    Q(:,iter+1) = Q(:,iter+1)./length(rt);
end

% calculate the return of the final rollout
Return(iter+1) = Q(1,iter+1);

% plot the return over the rollouts
figure(1);
plot(Return);
ylabel('return');
xlabel('rollouts');

disp(['Final Return ', num2str(Return(end))]);
