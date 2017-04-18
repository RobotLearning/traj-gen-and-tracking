%% Analyze ball data and train an adaptive (learning) filter
% that better predicts as we add more and more trials
%
% 1. Load all of the data and separate into trials
%    Check that the trials look sensible
% FOR each trial
% 2. Smoothen the data with a (Extended) Kalman Smoother
%    using current parameters
% 3. Data consists of bouncing balls hence
%    update flight parameters (incl. spin) and bounce params.
% 4. Check that the predictions are decreasing!

clc; clear; close all;
loadTennisTableValues;
file = '~/Dropbox/data/balls.txt';
M = dlmread(file);
% ball data
B = M(:,1:10); % cam id and status cols

[b1,b3,~] = prune_ball_data(B);
[balls1,balls3] = get_trials(b1,b3); %3d array

% plot one to see if it looks sensible
num_trials = length(balls1);
trial = 3;
plot_trial(balls1{trial}(:,2:4),balls3{trial}(:,2:4));

% Estimate parameters of a spin-model
spin.flag = false;
spin.Clift = Clift;

% estimate parameters before bounce
p = estimate_flight_params(balls3,spin);

% initialize EKF with a spin-model
% initialize ball with zero spin
% w0 = zeros(3,1);
% spin.flag = true;
% spin.Clift = Clift;
% spin.est = w0;
% 
% [ekfFilter,ball_func] = initFilterEKF(spin);

