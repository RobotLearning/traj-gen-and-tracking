%% Calculate probability of hitting 
% for a particular ball and racket trajectory and orientation
clc; clear; close all;

%% Compute the probability of hitting for a fixed racket trajectory 

% fix a racket trajectory
% fix a ball distribution (use kalman filter)
% see if simulations agree with the theoretical prob(hitting) = int(O to T)
% hitting time distribution

d = 2;
n = 100000;
T = 1.0;
dt = 0.02;
t = dt:dt:T;
N = length(t);
% generate a 2D ball trajectory
mu_b0 = [2; 2];
s_b0 = [0.01, 0.0; 0.0, 0.01];
mu_v0 = [-1; -1];
s_v0 = [0.1, 0.0; 0.0, 0.1];
% generate n sample trajectories
v0 = repmat(mu_v0,1,n) + chol(s_v0)*randn(2,n);
b0 = repmat(mu_b0,1,n) + chol(s_b0)*randn(2,n);
b = @(t) mu_b0 + mu_v0*t;
ball = zeros(d,N,n);
for i = 1:n
    ball(:,:,i) = repmat(b0(:,i),1,N) + v0(:,i)*t;
end

% generate a fixed racket trajectory
r1 = 0.1*rand; r2 = 0.1*rand;
w1 = 1 + 0.1*rand; w2 = 1 + 0.1*rand;

x = @(t) r1 + w1*t;
y = @(t) r2 + w2*t;
r = @(t) [x(t);y(t)];

xt = x(t);
yt = y(t);
rt = [x(t);y(t)];

% calculate the percentage of intersecting trajectories
% here we imagine a racket in the y-direction
disp('Calculating by direct simulation');
radius = 0.01;
dist = zeros(1,n);
tol = radius; %1e-2;
hit = 0;
for i = 1:n
    diff = rt - ball(:,:,i);
    % first get the index where they intersect
    %[M,I] = min(abs(diff'));
    % interpolate to get better accuracy
    dist = interp1(diff(1,:),diff(2,:),0,'spline');
    % look at the distance in the y-direction
    %dist = diff(2,I(2));
    %dist = min(diag(diff'*diff));
    if abs(dist) < tol
        hit = hit + 1;
    end
end
probHit = hit/n
        
% compute the probability of hitting
% rejection sampling here
%v0 = repmat(mu_v0,1,M) + chol(s_v0)*randn(2,M);
%b0 = repmat(mu_b0,1,M) + chol(s_b0)*randn(2,M);
disp('Rejection sampling')
Tcand = (b0(1,:) - r1)./(w1 - v0(1,:));
Thit = Tcand(abs(b0(2,:)-r2 + (v0(2,:)-w2).*Tcand)<radius);
%percHitTillT = sum(Thit < T)/length(Thit)
probHit = length(Thit)/n
probHitTillT = sum(Thit < T)/n
histogram(Thit)


disp('Integration of Gaussians')

mu1 = mu_b0(1) - r1;
mu2 = mu_b0(2) - r2;
nu1 = w1 - mu_v0(1);
nu2 = w2 - mu_v0(2);
s1 = s_b0(1,1);
s2 = s_b0(2,2);
g1 = s_v0(1,1);
g2 = s_v0(2,2);
tau1 = @(x,y) (-nu1.*(tol+x) - mu1.*y) ./ sqrt(g1.*(tol+x).^2 + s1.*y.^2);
tau2 = @(x,y) (nu1.*(tol-x) - mu1.*y) ./ sqrt(g1.*(tol-x).^2 + s1.*y.^2);
tauMin = @(x,y) min(tau1(x,y),tau2(x,y));
tauMax = @(x,y) max(tau1(x,y),tau2(x,y));
gaussSpace1 = @(x) exp(-(x-mu2).^2 ./ (2*s2)) ./ sqrt(2*pi*s2);
gaussSpace2 = @(y) exp(-(y-nu2).^2 ./ (2*g2)) ./ sqrt(2*pi*g2);
gaussTime = @(t) exp(-t.^2 ./ 2) ./ sqrt(2*pi);
fun = @(t,x,y) gaussTime(t).*gaussSpace1(x).*gaussSpace2(y);
probHit = integral3(@(X,Y,T)arrayfun(fun,T,X,Y),-Inf,Inf,-Inf,Inf,...
                    tauMin,tauMax)

% MU = @(t) b(t) - r(t);
% SIGMA = @(t) s_b0 + (s_v0 .* t^2);
% fun = @(t,x,y) exp(-(1/2)*(([x;y]-MU(t))'*(SIGMA(t)\([x;y]-MU(t))))) ./ ...
%        (((2*pi)^(d/2))*(det(SIGMA(t))^(1/2)));
% fun = @(t,y) fun(t,0,y);
% xmin = [0;-tol];
% xmax = [T;tol];
% tic;
% probHitTillT = integral2(@(T,Y)arrayfun(fun,T,Y),...
%           xmin(1),xmax(1),xmin(2),xmax(2));
% probHit = integral2(@(T,Y)arrayfun(fun,T,Y),...
%           0,Inf,-tol,tol);
% scale = integral2(@(T,Y)arrayfun(fun,T,Y),...
%            0,Inf,-Inf,Inf);
% %scale = 0.5;
% probHit = probHit / scale
% probHitTillT = probHitTillT / scale

% plot some samples
% b1 = squeeze(ball(1,:,:));
% b2 = squeeze(ball(2,:,:));
% plot(b1,b2)
% hold on;
% plot(xt,yt,'r--','LineWidth',2)
% hold off