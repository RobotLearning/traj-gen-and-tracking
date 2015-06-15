%% Testing total least squares


num_data = 500;
num_param = 10;
% design matrix
X = rand(num_data,num_param);
% corrupt design matrix
std_x = 0.1;
Xmeas = X + std_x * randn(size(X));
% hope that its rank is high
beta = rand(num_param,1);
% standard deviation of noise
std_y = 0.1;
y = X*beta + std_y * randn(num_data,1);

beta_ls = Xmeas\y;
beta_tls = tls(Xmeas,y,0.5);
disp('beta:');
disp(beta);
disp('beta_ls:');
disp(beta_ls);
disp('beta_tls:');
disp(beta_tls);
disp('norm of lsq error:');
disp(norm(beta-beta_ls,2));
disp('norm of tls error:');
disp(norm(beta-beta_tls,2));
