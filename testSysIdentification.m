%% Test system identification for ILC

% TODO: increase N,dimx,dimu and decrease trials. How does lsqlin compare?

clc; clear; close all

%% create A,B matrices and generate F
dimx = 1;
dimu = 1;

A = rand(dimx);
B = rand(dimx,dimu);
N = 20;
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

%% try to estimate Fact

% generate some trials
trials = 50;
U = rand(N*dimu,trials);
% get errors
E = Fact * U;

% estimate the F2 matrix
Fhat1 = (U' \ E')';
Ukron = kron(U,eye(N))';
FhatVec = Ukron \ E(:);
Fhat2 = reshape(FhatVec,N*dimx,N*dimu);

% estimation that preserves lower block diagonal structure of F
% get the indices that are zero
indZero = find(F == 0);
indNonZero = find(F ~= 0);
% drop those columns from Ukron
UkronDrop = Ukron(:,indNonZero);
FhatVecDrop = UkronDrop \ E(:);
FhatVec2 = zeros(N*N*dimx*dimu,1);
FhatVec2(indNonZero) = FhatVecDrop;
Fhat3 = reshape(FhatVec2,N*dimx,N*dimu);

% start x0 = F
Fvec = F(:); 
FvecDrop = Fvec(indNonZero);
FhatVecDrop2 = lsqlin(UkronDrop,E(:),[],[],[],[],[],[],FvecDrop);
FhatVec3 = zeros(N*N*dimx*dimu,1);
FhatVec3(indNonZero) = FhatVecDrop2;
Fhat4 = reshape(FhatVec3,N*dimx,N*dimu);

err = norm(Fhat1-Fact,'fro')

