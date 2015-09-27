%% Testing total least squares

clc; clear; close all;

%% Simple Test for TLS
% error gets worse as standard deviation on y increases

% m = 500; % number of data points
% n = 10; % number of variables
% % design matrix
% Z = rand(m,n);
% assert(rank(Z) == n, 'latent data matrix not full rank!');
% % corrupt design matrix
% std_x = 0.1;
% X = Z + std_x * randn(size(Z));
% % hope that its rank is high
% assert(rank(X) == n, 'corrupted matrix not full rank!');
% beta = rand(n,1);
% % standard deviation of noise
% std_y = 0.1;
% y = Z*beta + std_y * randn(m,1);
% 
% beta_ls = pinv(X,0.05)*y;
% beta_tls = tls(X,y,0.05);
% % disp('beta:');
% % disp(beta);
% % disp('beta_ls:');
% % disp(beta_ls);
% % disp('beta_tls:');
% % disp(beta_tls);
% disp('norm of lsq error:');
% disp(norm(beta-beta_ls,2));
% disp('norm of tls error:');
% disp(norm(beta-beta_tls,2));

%% create A,B matrices and generate F
dimx = 4;
dimu = 4;
% choose iteration number 
K = 10;

A = rand(dimx);
B = rand(dimx,dimu);
A2 = A;
B2 = B;
N = 20;
F = zeros(N*dimx,N*dimu);
Z = zeros(N*dimx,N*dimu,K);

% generate F
for i = 1:N
    for j = 1:i
        vec_x = (i-1)*dimx + 1:i*dimx;
        vec_u = (j-1)*dimu + 1:j*dimu;
        mat = B;
        for k = j+1:i
            mat = A * mat;
        end
        F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
    end
end

% generate multiple Z
for k = 1:K-1
    A2 = A2 + 0.001 * randn(dimx);
    B2 = B2 + 0.001 * randn(dimx,dimu);
    for i = 1:N
        for j = 1:i
            vec_x = (i-1)*dimx + 1:i*dimx;
            vec_u = (j-1)*dimu + 1:j*dimu;
            mat2 = B2;
            for m = j+1:i
                mat2 = A2 * mat2;
            end
            Z(vec_x,vec_u,k) = mat2;
        end
    end
end

% give info on matrices
% sigma_F = svd(F);
% sigma_Z = svd(Z);
% fprintf('Mismatch norm: %f\n', norm(F-Z,2));
% fprintf('sigma_max of F: %f\n', sigma_F(1));
% fprintf('sigma_max of Z: %f\n', sigma_Z(1));
% fprintf('sigma_min of F: %f\n', sigma_F(end));
% fprintf('sigma_min of Z: %f\n', sigma_Z(end));
% fprintf('Condition of F: %f\n', cond(F));
% fprintf('Condition of Z: %f\n', cond(Z));

% choose disturbance
d = 0.5 * rand(N*dimx,1);
fprintf('norm of disturbance: %f\n', norm(d,2));
u_ls = zeros(N*dimu,K);
e_ls = repmat(d,1,K);
e_tls = repmat(d,1,K);
u_tls = zeros(N*dimu,K);

for k = 1:K-1
    u_ls(:,k+1) = u_ls(:,k) - F \ e_ls(:,k);
    e_ls(:,k+1) = Z(:,:,k) * u_ls(:,k+1) + d;
    fprintf('norm of lsq error:%f\n',norm(e_ls(:,k+1),2));
    
    u_tls(:,k+1) = u_tls(:,k) - tls(F,e_tls(:,k),0.05);
    e_tls(:,k+1) = Z(:,:,k) * u_tls(:,k+1) + d;
    fprintf('norm of tls error:%f\n',norm(e_tls(:,k+1),2));
end
    
%disp('norm of minimum possible error:');
%e_min = (eye(size(Z,1)) - Z * pinv(Z)) * d;
%disp(norm(e_min,2));