% This function assumes you already run step0_plot_save_all_refs.m such
% that the file ref_reaching_all.mat was created and stored.
% 
% This script will load ref_reaching_all.mat and run ILC for each of the
% reaching points, independently.
% The number of ILC iterations are given by the parameter nMax.
%
% For each of the reaching points a corresponding file with the ILC
% solution will be saved as
%
%  ilcsol_rA_30.mat
%  ilcsol_rB_30.mat
%  ilcsol_rC_30.mat
%  ilcsol_rD_30.mat
%  ilcsol_rE_30.mat



%% initializing...
clear; clc; close all; dbstop if error;
add_required_folders();

%data_folder = fileparts(pwd);
data_folder = pwd;
data_folder = [data_folder '/common_data/'];

% default view angle [x y z] of the plot (this is just eye-candy)
dfv = [0.25 -1  0.25];
dfv = [0.0 -1   0.0];

% Create a "bioRob-based" robot here
[bioRob, cnf] = mdl_bioRob();


%% main code

% load ref_reaching_all.mat <= this file was previously created with the 
% script step0_plot_save_all_refs.m
load([data_folder 'ref_reaching_all.mat'])
nRefs = numel(refc);

% run this for loop for each of the reaching positions
for k= 1 : nRefs 
    
    prm.plot.refAnimation=1;
    
    fprintf('\n\n\n');
    fprintf('k: %g out of %g\n',k, nRefs);
    close all
    ref = refc{k};
    h = prepare_figures(dfv, ref, [], []);

    % animate the ideal tracking of the reference with zero errror
    if prm.plot.refAnimation
        figure(h.anim);
        title('Ideal tracking of the reference');
        animate_bioRob(bioRob, ref.q, 25);
    end

    % load controller parameters.
    % Without a second input argument, the feedforward inputs of the first
    % iteration will be set to zero.
    % Those feedforward inputs are in prm.ilc.u.
    prm = load_bioRob_controller(ref);
    
    % change to zero if you dont want to see animation after each iteration
    prm.plot.animation = 0; 
    
    % main learning loop    
    nMax   = 30; % number of ILC iterations
    [iter, prm]  = ilc_main_learning_loop([], nMax, prm, bioRob, h, ref);
   
    
    if 1 % enable saving
        refname = ref.name;
        filename = ['ilcsol_' refname '_' num2str( size(iter,2), 2 ) '.mat'];
        fprintf('Saving ILC solution with name: %s\n', filename)
        save([data_folder filename], 'iter', 'prm');
    end
   
end











