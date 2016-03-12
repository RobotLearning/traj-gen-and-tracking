%% Testing monotonic learning 

clc; clear; close all;

% 1d nonlinear system

% 
% eps = 0.2;
% dyn = @(x,u) x^2 + eps*x - u;
% nomDyn = @(x,u) x^2 - u;
% 
% a = 2;
% unom = @(x) x^2 - a;
% x0 = 0;
% fnc = @(t,x) dyn(x,unom(x));
% [T,X] = ode45(fnc,[0 1], x0);
% 
% N = 100;
% t = linspace(0,1,N);
% Xt = interp1(T,X,t);
% plot(t,a*t,'--',t,Xt,'-r');

%% Estimate linear system with least squares

% assumption is that the trials are exactly the same
% for linear time varying system F
% and random inputs

%{
m = 2;
n = 4;
N = 4;
K = 8; % number of trials has to be N*m
B = rand(n,m,N);
A = rand(n,n,N);
F = liftDyn(A,B);
u = rand(N*m,K);
e = F*u;

% estimate F
for i = 1:K
    U((N*(i-1))+1:(N*(i-1))+N,1:m*N*N) = kron(eye(N),u(:,i)');
    E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
end
D = duplicateVech(N,m);
M = U * D';
Est = pinv(M) * E;
EstAddZero = Est' * D;

for i = 1:N
    F_est((i-1)*n+1:i*n,:) = EstAddZero(:,N*m*(i-1)+1:N*m*i);
end

errNorm = norm(F-F_est)
%}


%% Apply plant inversion ILC

%{
n = 1; % dim_x
m = 1; % dim_u
T = 1.0; % final time
N = 50;
dt = T/N; % discretization
A = rand(n,n,N);
B = rand(n,m,N);
% TODO: make sure they are correlated not purely random
A_act = A + 0.001*rand(n,n,N);
B_act = A + 0.001*rand(n,m,N);

[Ad,Bd] = discretizeDyn(A,B,dt);
[Ad_act,Bd_act] = discretizeDyn(A_act,B_act,dt);
% the unknown system
Z = liftDyn(Ad_act,Bd_act);
% the nominal system
F = liftDyn(Ad,Bd); %0.80 * F_act; %liftDyn(Ad,Bd);

% theoretical computations
% try
%     s_nom = svd(F,'econ');
%     s_diff = svd(Z-F,'econ');
%     s_max_nom = max(s_nom);
%     s_min_nom = min(s_nom);
%     s_max_diff = max(s_diff);
%     mu = (1 + sqrt(5))/2;
%     BOUND = s_min_nom^2 / (mu * s_max_nom + mu * s_min_nom + s_min_nom);
%     spectral_radius = abs(max(eig(eye(size(Z,2)) - pinv(F)*Z)));
%     assert(s_max_diff < BOUND,'no monotonic convergence!');
%     assert(spectral_radius < 1, 'no asymptotic convergence!');
% catch ME
%     disp(ME.message);
% end
    

% observe errors
numTrials = 100; % trials 
err = zeros(N,numTrials);
x0 = zeros(n,1);
us = zeros(N,numTrials+1);
err_norm = zeros(1,numTrials);

% reference traj
a = 2;
t = linspace(0.02,1,N);
ref = a*t(:);

alpha = 0.1;
% we assume no noise for now
for i = 1:numTrials
    % get error
    err(:,i) = Z * us(:,i) - ref;
    us(:,i+1) = us(:,i) - pinv(F,0.05) * err(:,i);
    err_norm(i) = norm(err(:,i),2);
end

scatter(1:numTrials,err_norm);
%}

%% Apply adaptive ilc

n = 1; % dim_x
m = 1; % dim_u
T = 1.0; % final time
N = 20;
dt = T/N; % discretization
A = rand(n,n,N);
B = rand(n,m,N);
% TODO: make sure they are correlated not purely random
A_act = A + 0.001*rand(n,n,N);
B_act = A + 0.001*rand(n,m,N);

[Ad,Bd] = discretizeDyn(A,B,dt);
[Ad_act,Bd_act] = discretizeDyn(A_act,B_act,dt);
% the unknown system
Z = liftDyn(Ad_act,Bd_act);
% the nominal system
F = liftDyn(Ad,Bd); %0.80 * F_act; %liftDyn(Ad,Bd);    

% observe errors
numTrials = 10; % trials 
err = zeros(N,numTrials);
x0 = zeros(n,1);
us = zeros(N,numTrials+1);
err_norm = zeros(1,numTrials);

% reference traj
a = 2;
t = linspace(0.02,1,N);
ref = a*t(:);

% initialize Gamma matrix and a, b 
var_nom = 1;
a = 2; 
b = var_nom;
Gamma = 100.0 * eye(n); % prior precision matrix, akin to saying we have 100 data points

% we assume no noise for now
for i = 1:numTrials
    % get error
    err(:,i) = Z * us(:,i) - ref;
    % estimate system
    F = estimatePlantBayes(F,G,a,b,us(:,1:i),err(:,1:i),N);
    us(:,i+1) = us(:,i) - pinv(F_est,0.05) * err(:,i);
    err_norm(i) = norm(err(:,i),2);
end

scatter(1:numTrials,err_norm);