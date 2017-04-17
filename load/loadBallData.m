
function [t,B] = loadBallData(dataSet)

dataSet = 1;
file = ['~/Dropbox/data/realBallsimRobotData',int2str(dataSet)];
M = dlmread([file,'.txt']);
% ball data
B = M(1:2:end,:);
% robot joint data
N_DOFS = 7;
Q = M(2:2:end,1:2*N_DOFS);
scale = 0.001; % recorded in milliseconds
t = scale * B(:,11);

end