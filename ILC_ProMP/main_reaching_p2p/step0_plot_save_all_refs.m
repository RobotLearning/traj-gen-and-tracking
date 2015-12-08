%
%
%
% This file is used to load a pseudo-bioRob and find inverse kinematics (IK) for
% different reaching positions. Those positions are named rA, rB, rC rD rE.
% The parameters of the reaching positions are pre-computed and stored in
% the file create_reference.m
% 
% This script will also compute the inverse dynamics (ID) solution as a ground
% truth for ideal torques to reach the rA...rE positions. Those torques can
% be used in the future for comparison of the proposed method.
%



clear; clc; close all; dbstop if error;
add_required_folders();

%data_folder = fileparts(pwd);
data_folder = pwd;
data_folder = [data_folder '/common_data/'];

% default view angle [x y z] of the plot (this is just eye-candy)
dfv = [0.1 -1  0.1];

% Create a "pseudo-bioRob" robot here
[bioRob, cnf] = mdl_bioRob();

prm.plot.refAnimation=1;


%% main code

% ===========================
% Create trajectories without plotting them.
% Find the IK and ID for different reaching positions.
% ==========================
if 0
    refnameS = [{'rA', 'rB', 'rC', 'rD', 'rE'}];
    for m = 1:numel(refnameS)
        refname       = refnameS{m};
        ref_file_name = ['ref_' refname '.mat'];
        refc{m}       = create_reference(bioRob, 0.02, refname);
        ref           = refc{m};
        save([data_folder ref_file_name], 'ref');
    end
    % store the reaching positions as a mat file.
    save([data_folder 'ref_reaching_all.mat'], 'refc');
 end

% ===========================
% Load and plot all reference trajectories
% ==========================
load([data_folder 'ref_reaching_all.mat']);
h = plot_all_references(dfv, bioRob, refc);


%           sspng([data_folder], 'reach')














