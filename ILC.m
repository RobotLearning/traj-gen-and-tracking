function [err,u_app] = ILC(t,x_nom,u_trj,fun,STR,num_iters)

%%%%%%%%%%%%%%%%%% ILC CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs :
% 
% t - time vector
% x_nom - nominal trajectory
% u_trj - nominal control input
% fun - function handle for generating state errors
% STR - structure containing necessary learning parameters/matrices
% num_iter - number of iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% RELEASE STRUCTURE %%%%%%%%%%%%%%%%%%%%%%
dim_x = STR.dim_x;
dim_u = STR.dim_u;
w = STR.w;
lsq_cost = STR.lsq_cost;
PAR = STR.PAR;
COV = STR.COV;
CON = STR.CON;
% Lifted-domain matrices for Kalman filter and convex optimization
F = STR.F;
G = STR.G;
H = STR.H;
% scaling matrix
S = STR.S;
% constraints
umin = STR.umin; 
umax = STR.umax;
% covariances
Omega = STR.COV.Omega; % process noise covariance
M = STR.COV.M; % covariance matrices to be passed for calc. error
P0 = STR.P0;
d = STR.d0; % estimated error
P = P0; % variance of the error estimate

%%%%%%%%%%%%% INITIALIZE VARIABLES %%%%%%%%%%%%%%
iter = 1;
N = length(t);
Nu = length(t) - 1;
h = t(2) - t(1); % assuming uniform discretization
% holds the trajectory deviation norms
diff = zeros(1,num_iters); 
% deviation of control input
u = zeros(dim_u*Nu,1); 
% store applied inputs here
u_app = zeros(dim_u, N, num_iters); 
% store errors (traj. deviation) here
err = Inf(dim_x*N,num_iters); 

%%%%%%%%%%%%% START ITERATION %%%%%%%%%%%%%%%%%%%
while iter <= num_iters
    
    fprintf('Iteration no %d \n', iter); 
    % complete u to match the dimensions of x
    u_iter = [u; u(end-dim_u+1:end)]; 
    % regroup into suitable matrix for simulation
    u_in = reshape(u_iter,dim_u,N);
    % save for later use
    u_app(:,:,iter) = u_trj + u_in;
    % get error (observed trajectory deviation)
    [x_iter, y_dev] = fun(t, u_trj, u_in, x_nom, PAR, COV, CON);
    % get frobenius-norm error
    diff(iter) = lsq_cost(x_iter,x_nom,w);
    err(:,iter) = y_dev;
    
    %%%%%%%%%%%% STATE ERROR ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%
    % estimate error d with kalman filter
    % calculate Kalman gain K and estimate of error d
    % Omega can decay to avoid adapting to outliers and nonrepetitive noise
    % Omega = Omega * 0.95;
    P = P + Omega;  
    Theta = G * P * G' + M;
    K = P * G' * inv(Theta); %#ok
    P = P - K * G * P;
    d = d + K * (y_dev - G*d - (G*F + H)*u);

    %%%%%%%%%%%%% FORM CONSTRAINT AND PENALTY MATRICES %%%%%%%%
    % arrange in format Lu <= q
    % call the particular constraints-generating code
    if ~isfield(STR, 'L'), 
        [L,q] = cnstr(STR,d); 
    else L = STR.L; q = STR.q;
    end
    % input deviation penalty matrix D
    D0 = eye(length(u)); 
    D1 = (diag(ones(1,dim_u*(Nu-1)),dim_u) - eye(dim_u*Nu))/h;
    D1 = D1(1:end-dim_u,:); % D1 is (Nu-1)*nu x Nu*nu dimensional
    D2 = (diag(ones(1,dim_u*(Nu-2)),2*dim_u) - ... 
         2*diag(ones(1,dim_u*(Nu-1)),dim_u) + eye(dim_u*Nu))/(h^2);
    D2 = D2(1:end-2*dim_u,:); % D2 is (Nu-2)*nu x Nu*nu dimensional
    % penalty scale
    a0 = 5e-5; a1 = 5e-5; a2 = 5e-5;
    % slack for input u
    eps_u = 1e-6;
    
    %%%%%%%%%%%%%%%%%% SOLVE OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%
    if iter < num_iters
    % solve with quadprog
    options = optimset('Display', 'iter', 'Algorithm', 'interior-point-convex');
    u = quadprog(2*(F'*S')*S*F + 2*a2*(D2'*D2), 2*F'*S'*d, L, q, [], [], umin(:), umax(:), [], options);
    
    % call convex optimization with CVX
    % find u which minimizes state error x_dev
%     cvx_begin
%         % number of variables to optimize is N 
%         l = 2; % 2,1,or Inf norms can be considered
%         variable u(u_dim*N) % control input u
%         minimize(norm(S*(F*u + d),l) + a0*norm(D0*u,l))
%         subject to 
%             umin(:)+eps_u <= u <= umax(:)-eps_u %#ok
%             L*u <= q %#ok
%     cvx_end
    end
    
    % plot iteration results
    ilc_results(t,x_iter,u_app,iter);
    iter = iter + 1; 
end

%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot final x trajectory
fprintf('Total of %d iterations \n', iter-1);
figure; plot(t,x_iter);
title('The final trajectories');

% Plot convergence
figure; plot(1:num_iters,diff);
title('Trajectory deviation 2-norm error');