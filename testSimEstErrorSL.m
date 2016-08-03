% Test Estimation Error in SL simulation (using KF and nonlinear ball model)

% load the file
loadTennisTableValues;
time = '3-8-16_14:46:24';
file = ['../robolab/barrett/saveData/sim_data_',time];
M = dlmread([file,'.txt']);
% ball data
B = M(:,1:1+6+6);
t = B(:,1)/1000;
ballSim = B(:,2:7);
ballPred = B(:,8:end);

% find the times at which there is a jump in Y - reset times
thresh = 0.8;
diffBalls = diff(ballSim(:,1:3));
ballDist = diag(diffBalls*diffBalls');
idx = find(ballDist > thresh);
resetIdx = [1;idx];

% get the index at which robot hits the ball
idxHit = find(ballSim(1:end-1,5) .* ballSim(2:end,5) < 0 & ballSim(1:end-1,5) > 0);

% plot one case
trial = 1;
idxPlot = resetIdx(trial):idxHit(trial);

% plot the case where there was a hit
figure;
scatter3(ballSim(idxPlot,1),ballSim(idxPlot,2),ballSim(idxPlot,3),'r');
hold on;
scatter3(ballPred(idxPlot,1),ballPred(idxPlot,2),ballPred(idxPlot,3),'b');
title('Predicted ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
% s1.MarkerEdgeColor = s1.CData;
% s2.MarkerEdgeColor = s2.CData;
legend('sim ball','ball pred');
hold off;

% plot the estimation error
err = ballSim(idxPlot,:) - ballPred(idxPlot,:);
err_norm = diag(err*err');
figure;
min_obs_to_start_filter = 100;
plot(err_norm(min_obs_to_start_filter+1:end));
ylabel('Estimation error norm');
xlabel('Number of observations');