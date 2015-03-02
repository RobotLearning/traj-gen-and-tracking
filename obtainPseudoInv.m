%% Construct pseudoinverse of matrix loaded from a file

% SL dumps the matrix and loads its pseudoinverse

folder = '../robolab/barrett/prefs/';
file = [folder,'F.txt'];
F = dlmread(file);
dof = 7;

len = size(F,2)/dof; %length of the trajectory
%Finv = pinv(F);

Q = diag([ones(1,4), 0.1*ones(1,3), 0.0*ones(1,4), 0.0*ones(1,3)]);
R = 0.01 * diag([ones(1,4),zeros(1,3)]);

% deal Q matrix to Ql
Ql = cell(1,len);
Rl = cell(1,len);

[Ql{:}] = deal(Q);
[Rl{:}] = deal(R);
Ql = blkdiag(Ql{:});
Rl = blkdiag(Rl{:});

Finv = (F' * Ql * F + Rl) \ (F' * Ql);
max(max(Finv))

dlmwrite([folder,'Finv.txt'], Finv, 'delimiter', '\t');