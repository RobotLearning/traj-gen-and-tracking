%% Script that plots the trajectories in SL

cd('/home/okankoc/');
% extract the desired positions
folder = 'robolab/barrett/saveData/';
fileName = 'joint_des_pos.txt';

% load the file
file = [folder,fileName];
M = dlmread(file);
h = 0.002;
t = h * (1:size(M,1));

for i = 1:7
    names{i} = ['joint_', int2str(i)];
    err{i} = ['err_j_', int2str(i)];
end

% Compare with actual trajectories
fileNameAct = 'joint_pos.txt';
fileAct = [folder,fileNameAct];
Mact = dlmread(fileAct);

figure;
plot(t,M,'--',t,Mact,'-');
legend(names);

figure;
plot(t,M-Mact);
legend(err);

