%% Script that plots the trajectories in SL

clc; clear; close all;
cd('/home/okankoc/robolab/barrett/saveData/');
% extract the desired positions
%folder = 'robolab/barrett/saveData/';
fileName = 'joint_des_pos.txt';

% load the file
%file = [folder,fileName];
M = dlmread(fileName);
h = 0.002;
t = h * (1:size(M,1));
qdes = M(:,2:end);

for i = 1:7
    names{i} = ['joint_', int2str(i)];
    err{i} = ['err_j_', int2str(i)];
end

% Compare with actual trajectories
fileNameAct = 'joint_act_pos.txt';
%fileAct = [folder,fileNameAct];
Mact = dlmread(fileNameAct);
qact = Mact(:,2:end);

figure(1);
plot(t,qact);
legend(names);

figure(2);
plot(t,qdes,'--',t,qact,'-');
%legend(names);

figure(3);
plot(t,qact-qdes);
legend(err);

