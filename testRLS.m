%% Test recursive least squares for ILC

clc; clear; close all
% create A,B matrices
dimx = 4;
dimu = 2;

A = rand(dimx);
B = rand(dimx,dimu);
N = 100;
F = zeros(N*dimx,N*dimu);

Fact = zeros(N*dimx,N*dimu);
A2 = A + 0.01 * randn(dimx);
B2 = B + 0.01 * randn(dimx,dimu);

% generate F
for i = 1:N
    for j = 1:i
        vec_x = (i-1)*dimx + 1:i*dimx;
        vec_u = (j-1)*dimu + 1:j*dimu;
        mat = B;
        mat2 = B2;
        for k = j+1:i
            mat = A * mat;
            mat2 = A2 * mat2;
        end
        F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
        Fact(vec_x,vec_u) = mat2;
    end
end

% try to estimate F2
% generate some trials
trials = 10;
U = rand(N*dimu,trials);
% get errors
E = Fact * U;

% estimate the F2 matrix
S = cell(1,trials);
[S{:}] = deal(U);
Ubar = blkdiag(S{:});
F_hat_vec = Ubar \ E(:);
F_hat = reshape(F_hat_vec,N*dimx,N*dimu);

err = norm(F_hat,F2)

