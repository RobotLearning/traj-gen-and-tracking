% This script will run ILC with initial guesses provided by the previous
% step.
% It should be obvious by the plot of the RMS error that ProMPs with a
% large number of ILC trajectories will provide better prediction.
% In fact, when you use the ProMP created from the ILC training data of
% rA+rB+rC+rD the initial tracking to reach rE is already close to the
% baseline error.
% 


clear; clc; close all; dbstop if error;

add_required_folders();

data_folder = fileparts(pwd);
data_folder = [data_folder '/common_data/'];
    

% default view angle [x y z]
dfv = [0.1 -1  0.1];


% Create a "bioRob-based" robot here
[bioRob, cnf] = mdl_bioRob();


%% main code

refname = 'ilc_conditioned_guess_all_reach.mat';
load([data_folder refname]);

nCases = numel(iguess);
cl = distinguishable_colors(nCases);

h =[];
for k = 1 : nCases
    
    prm.plot.refAnimation=0;
    
    fprintf('\n\n\n');
    fprintf('k: %g out of %g\n',k, nCases);
    %close all
    
    ref = iguess{k}.ref;
    h = prepare_figures(dfv, ref, [], h);
    h.colors = cl(k,:);

    % animate the ideal tracking of the reference with zero errror
    if prm.plot.refAnimation
        figure(h.anim);
        title('Ideal tracking of the reference');
        animate_bioRob(bioRob, ref.q, 10);
    end

    % load controller parameters
    if isempty(iguess{k}.training_data) % this means that no promp was used
                                        % but the initial guess came
                                        % directly from some other case
        prm = load_bioRob_controller(ref, iguess{k}.init_guess.u);
    else
        prm = load_bioRob_controller(ref, iguess{k}.init_guess.u_conditioned);
    end

    % change to zero if you dont want to see animation after each iteration
    prm.plot.animation = 1; 

    % main learning loop    
    % keep this 1 as default
    nMax   = 5;
    prm.comment = iguess{k}.comment;
    [iter, prm]  = ilc_main_learning_loop([], nMax, prm, bioRob, h, ref);
    
    final_ilc{k} = iguess{k};
    final_ilc{k}.final.iter = iter;
    final_ilc{k}.final.prm  = prm;   
end

if 0
    refname = ref.name;
    filename = ['ilc_final_' refname '_' num2str( size(iter,2), 2 ) '.mat'];    
    save(filename, 'final_ilc');
end















