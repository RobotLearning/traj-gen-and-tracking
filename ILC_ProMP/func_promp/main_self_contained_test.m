clear; close all; clc; dbstop if error;
try addpath('~/Dropbox/myMatlabFunctions/'); end
try addpath('localdata/Dropbox/myMatlabFunctions/'); end
addpath('func_process_data/');
addpath('func_plot/');
addpath('func_promp/');

strial{1} = struct('Color', 'b', 'LineStyle', '-', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','b', 'MarkerSize', 5);
strial{2} = struct('Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','r', 'MarkerSize', 5);
strial{3} = struct('Color', 'g', 'LineStyle', '-', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','g', 'MarkerSize', 5);
strial{4} = struct('Color', 'm', 'LineStyle', '-', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','m', 'MarkerSize', 5);

%% User inputs

% Loading data
% In both TYPE1 and TYPE2 data will return a system with 5 joints.
% q(1), q(2), q(3) represent the first 3 links of the arm
% q(4) and q(5) represent the X and Y position in Cartesian space

prm.data.q_qdot = 'q';
prm.data.u_ilc  = 'ilc';
prm.data.dt     = 0.01;

prm.nJointA=3;
prm.nJointB=3;

prm.plot.check_encoding = 0; 
prm.plot.trainData      = 0;
prm.plot.sameYscale     = 0;
prm.plot.nPlots = 5;

prm.obs.nTimes   =  20; % number of observations
prm.obs.stdev    =  0.1;
prm.obs.jointPos = [ 1 2 3 ]; % list here the joint positions  that are going to be observed
prm.obs.jointVel = [ 1 2 3 ];

observed_data_index = [ 1:5:100 ]; % as percentages


%% Code starts from here

data = [];

% This will load data set ilc_mid_10 and extract only the data_tst as a
% test data
[data, data3_aux]   = load_bioRobData_incremental('ilc_mid_20.mat',  data, prm);

if 0 % extract one of the iterations from ILC
    [data, data_tst]  = remove_test_data(data, 5);
else % use the reference itself as the test data
    data_tst = data3_aux.ref_resampled; 
end

data =[]; % get rid of all training related to the test data
[data, data1_aux]   = load_bioRobData_incremental('ilc_up_20.mat',    data, prm);
[data, data2_aux]   = load_bioRobData_incremental('ilc_down_20.mat',  data, prm);


[data]   = get_statistics(data);
data_ias = format_ias_data(data);


% plot individual joints
tmp.work_format= data;
tmp.taskA = 'ilc_bioRob';

time = data_tst(1).t;



h_error = figurew('error');
for trial = 1:1

    
    if trial==1
        prm.obs.jointPos = [ 1 2 3 ];
        prm.obs.jointVel = [ 1 2 3 ];
    end
    if trial==2
        prm.obs.jointPos = [ 1 2 3 ];
        prm.obs.jointVel = [ ];
    end
    if trial==3
        prm.obs.jointPos = [ ];
        prm.obs.jointVel = [ 1 2 3 ];
    end
    
    
    % first use of time is given here
    h = plot_nStates(time, tmp, prm);


    % create promp here
    pmp = regression_pmp(data_ias, 30, time);

    pmp.h = h;

    if  prm.plot.check_encoding  % check if trajectory was properly encoded by ProMP
        check_ProMP_encoding(pmp, data);
    end

    observed_data_index2 = uint32( observed_data_index/100*pmp.nTraj );
    if observed_data_index2(end)~=pmp.nTraj
        observed_data_index2(end+1)=pmp.nTraj;
    end

    obs = Observation(prm.obs.stdev, prm, data_tst); % CREATE object
    obs.measuredIndexes(observed_data_index2);

    % ================================
    % plot test data
    % ================================
    figure(pmp.h.qA);
    for j=1:prm.nJointA % position
        subplot(1, prm.nJointA,j);
        plot(time, obs.groundtr(j).q, struct('Color', 'k', 'LineStyle', '--', 'LineWidth', 1));
    end
    figure(pmp.h.qAdot);
    for j=1:prm.nJointA  % velocity
        subplot(1, prm.nJointA,j);
        plot(time, obs.groundtr(j).qdot, struct('Color', 'k', 'LineStyle', '--', 'LineWidth', 1));
    end
    figure(pmp.h.qB);
    for j=1:prm.nJointB
        subplot(1, prm.nJointA,j);
        plot(time, obs.groundtr(prm.nJointA+j).q, struct('Color', 'k', 'LineStyle', '--', 'LineWidth', 1));
    end

    [visited, joint_KF] = mainloop(pmp, prm, obs);

    figure(h_error);
    plot(time(visited.k), visited.mse, strial{trial});
    plot2svg('error')
    

end










































