%% Recursive least squares

n = 5;
N = 10;
X = rand(N,n);
beta = rand(n,1);
eps = randn(N,1);
y = X * beta + eps;

% estimate beta first
beta_est = X \ y;

% get a new data point
x_add = rand(1,n);
y_add = x_add * beta + randn(1);
X_new = [X; x_add];
beta_est_batch = X_new \ [y;y_add];

% apply recursive least squares
A = X'*X
alpha = 1/(1 + x_add*A*x_add');
C = eye(N) - alpha*A*B;
beta_est_rls = C * beta_est + delta;