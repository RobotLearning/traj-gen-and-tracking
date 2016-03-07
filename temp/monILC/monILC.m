%% Testing monotonic learning in 1d nonlinear system

% clc; clear; close all;
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

%% 2 by 2 linear system

clc; clear; close all;
n = 2;
N = 3;
K = 4; % number of trials
B = rand(n);
A = rand(n);
F = [B, zeros(n), zeros(n);
     A*B, B, zeros(n);
     zeros(n,2*n), B];
u = rand(N*n,K);
e = F*u;

% estimate F
for i = 1:K
    U((N*(i-1))+1:(N*(i-1))+N,1:n*N^2) = kron(eye(N),u(:,i)');
    E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
end
D = duplicateVech(N,n,n);
M = U * D';
Est = pinv(M) * E;
EstAddZero = Est' * D;
% FIXME: this part should be automated
F_est = [EstAddZero(:,1:N*n); EstAddZero(:,N*n+1:2*N*n)]
errNorm = norm(F-F_est)

%% Performance for a linear system

%{
clc; clear; close all;
n = 1; % dim_x
m = 1; % dim_u
T = 1.0; % final time
N = 50;
dt = T/N; % discretization
A = randn(n,n,N);
B = randn(n,m,N);
% TODO: make sure they are correlated not purely random
A_act = A - 0.01;
B_act = A - 0.01;

[Ad,Bd] = discretizeDyn(A,B,dt);
[Ad_act,Bd_act] = discretizeDyn(A_act,B_act,dt);
F = liftDyn(Ad,Bd);
F_act = liftDyn(Ad_act,Bd_act);

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

alpha = 0.1;
% we assume no noise for now
for i = 1:numTrials
    % get error
    err(:,i) = F_act * us(:,i) - ref;
    us(:,i+1) = us(:,i) + alpha* err(:,i);
    err_norm(i) = norm(err(:,i),2);
end

plot(err_norm);

%% Estimate system dynamics along traj

% produce duplication matrix
D = duplicateVech(N,n,m);

% assuming d is zero
Ubar = us(:,1:numTrials);
M = kron(Ubar,eye(N*m))' * D;
E = err;
vecE = E(:);
vechF = M \ vecE;
% duplicate to vecF
vecF = D * vechF;
% duplicate back to F
F_est = reshapeF(vecF,n,m,N);
%}