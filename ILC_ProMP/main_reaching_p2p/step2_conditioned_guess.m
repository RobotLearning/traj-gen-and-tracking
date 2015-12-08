%
% This file will create different ProMPs and condition each of them on the
% same desired trajectory given by the reaching point rE. The conditioned
% solution is a vector of joint torques that can be used to reach the position rE.
%
% There will be a total of 4 ProMPs.
% The first one uses only the trajectories obtained from the ILC when
% trying to track the trajectory to reach the position rA.
% The second ProMP uses the trajectories obtained from the ILC when
% trying to track the trajectory to reach the position rA and rB.
% The second ProMP uses the trajectories obtained from the ILC when
% trying to track the trajectory to reach the position rA, rB, rC.
% The second ProMP uses the trajectories obtained from the ILC when
% trying to track the trajectory to reach the position rA, rB, rC, rD.
%
% We expect that by adding more and more trajectories when learning the
% ProMP the prediction of the torques to reach rE will become better.
%
% Those ProMPs are then conditioned on the trajectory for rE, and the
% output is a vector of torque commands for each of the joint of the arm.
% This script uses only 21 equaly spaced points of the desired trajectory
% rE, rather than the whole trajectory. This is simply to make the code
% faster, but implementing the updates in batch will eliminate this
% limitation.
%
% Additionally, for the sake of comparison I also added a fifth solution
% where the solution is given by the torques to reach the point rD. This
% solution is used to evaluate a naive method of the type nearest-neighbor.
%
%
clear; clc; close all; dbstop if error;

extraFolder1 = add_required_folders_promp();

data_folder = fileparts(pwd);
data_folder = [data_folder '/common_data/'];
    
% default view angle [x y z] of the plot (this is just eye-candy)
dfv = [0.1 -1  0.1];

% Create a "bioRob-based" robot here
[bioRob, cnf] = mdl_bioRob();

%% start main code here

% load all ilc solutions that were learned independently in the script 
% step1_train_from_zero.m
isol = load_reaching_solutions(data_folder);

% desiredRef = 5 defines the reaching position rE.
desiredRef = 5; % 

% training data: rA
kk = 1;
iguess{kk}.training_data = 1;
iguess{kk}.comment = 'init_guess comes from promp trained with rA';
[iguess{kk}.init_guess, iguess{kk}.prompsol] = promp_conditioning(iguess{kk}.training_data,...
                      isol, extraFolder1, desiredRef, 'trainingData: rA');
iguess{kk}.ref = isol{desiredRef}.prm.ref;

% training data: rA+rB
kk = 2;
iguess{kk}.training_data = [1 2];
iguess{kk}.comment = 'init_guess comes from promp  trained with rA+rB';
[iguess{kk}.init_guess, iguess{kk}.prompsol] = promp_conditioning(iguess{kk}.training_data,...
                       isol, extraFolder1, desiredRef, 'trainingData: rAB');
iguess{kk}.ref = isol{desiredRef}.prm.ref;

% training data: rA+rB+rC
kk = 3;
iguess{kk}.training_data = [1 2 3];
iguess{kk}.comment = 'init_guess comes from promp trained with rA+rB+rC';
[iguess{kk}.init_guess, iguess{kk}.prompsol] = promp_conditioning(iguess{kk}.training_data,...
                       isol, extraFolder1, desiredRef, 'trainingData: rABC');
iguess{kk}.ref = isol{desiredRef}.prm.ref;

% training data: rA+rB+rC+rD
kk = 4;
iguess{kk}.training_data = [1 2 3 4];
iguess{kk}.comment = 'init_guess comes from promp trained with rA+rB+rC+rD';
[iguess{kk}.init_guess, iguess{kk}.prompsol] = promp_conditioning(iguess{kk}.training_data,...
                      isol, extraFolder1, desiredRef, 'trainingData: rABCD');
iguess{kk}.ref = isol{desiredRef}.prm.ref;


% training data: use the previous reaching posisiton rD as the initial
% guess to reach rE
kk = 5;
iguess{kk}.training_data = [];
iguess{kk}.comment    = 'use previous solution from rD as initial guess';
iguess{kk}.init_guess = isol{4}.prm.ref;
iguess{kk}.ref        = isol{desiredRef}.prm.ref;


if 1
    save([ data_folder  'ilc_conditioned_guess_all_reach.mat'], 'iguess');
end

if 0 % save plots of the conditioned promp
    sspng('plots_promps', []);
end







