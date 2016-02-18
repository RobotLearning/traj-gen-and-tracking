%% 3D table tennis with 2 Barrett WAMs

clc; clear; close all; dbstop if error;

% VHP method vs. optimal control approach

profile on;
initializeWAM;
% change base for wam2
wam2 = [];
% For a match run
% tt.match();
tt = TableTennis(wam,wam2);
tt.match();
profile viewer


