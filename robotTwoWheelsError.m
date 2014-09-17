function [traj, err] = robotTwoWheelsError(t, u, x, PAR, Omega, M, CON)

% calculates the swing-up error of the pendulum dynamics

dim = size(x,1);
% state trajectories for this iteration
traj = x;

par.wheel1.radius = 1.1 * PAR.wheel1.radius;
par.wheel2.radius = 0.9 * PAR.wheel2.radius;
par.length = PAR.length;

% simulate real trajectory with RK4
handle = @robotTwoWheelsKinematics;
traj = sim_RK4(t,traj,u,CON,par,handle);

% vectorize x_iter into N*dim dimensions
x_iter_vec = traj(:);
x_iter_vec = x_iter_vec(dim+1:end); % start from x(2)
% add lifted-domain error noise w with covariance Omega
x_vec_noisy = x_iter_vec + chol(Omega)*randn(length(x_iter_vec),1);
% add observation noise with covariance M
y = x_vec_noisy + chol(M)*randn(length(x_iter_vec),1);

% find deviation in observation
x_vec = x(:); 
x_vec = x_vec(dim+1:end);
err = y - x_vec;
% start from x(1)
%err = [zeros(dim,1); err];