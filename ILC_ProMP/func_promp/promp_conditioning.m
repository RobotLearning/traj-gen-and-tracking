function [ref, visited]   = promp_conditioning( trainingData, datac, path1, desRefNum, optionalName )

    %% initializing
    
%     extra_folders = ({'func_process_data/', 'func_plot/', 'func_promp/'} );
% 
%     % adding path so that promp works
%     for k=1:numel(extra_folders)
%         addpath([path1 '/' extra_folders{k}]);
%     end

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
    prm.plot.trainData      = 1;
    prm.plot.sameYscale     = 0;
    prm.plot.nPlots = 1;

    prm.obs.nTimes   =  20; % number of observations
    prm.obs.stdev    =  0.1;
    prm.obs.jointPos = [ 1 2 3 ]; % list here the joint positions  that are going to be observed
    prm.obs.jointVel = [ 1 2 3 ];

    observed_data_index = [ 1:5:100 ]; % as percentages


    %% Code starts from here
    
    % use this to extract the reference
    [data, data_aux] = format_data(datac{desRefNum}, prm.data.dt, [], prm.data.q_qdot, prm.data.u_ilc);
    data_tst = data_aux.ref_resampled;
    data = [];
    
    for k=1:numel(trainingData)
        data = format_data(datac{ trainingData(k) }, prm.data.dt, data, prm.data.q_qdot, prm.data.u_ilc);
    end

    data = get_statistics(data);
    data_ias = format_ias_data(data);


    % plot individual joints
    tmp.work_format= data;
    tmp.taskA = 'ilc_bioRob';

    time = data_tst(1).t;

    h_error = figurew('error');

    trial = 1;
    
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
    h = plot_nStates(time, tmp, prm, optionalName);


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
    legend({'training data', 'promp prior', 'ground truth'});
    
    figure(pmp.h.qAdot);
    for j=1:prm.nJointA  % velocity
        subplot(1, prm.nJointA,j);
        plot(time, obs.groundtr(j).qdot, struct('Color', 'k', 'LineStyle', '--', 'LineWidth', 1));
    end
    legend({'training data', 'promp prior', 'ground truth'});
    
    figure(pmp.h.qB);
    for j=1:prm.nJointB
        subplot(1, prm.nJointA,j);
        plot(time, obs.groundtr(prm.nJointA+j).q, struct('Color', 'k', 'LineStyle', '--', 'LineWidth', 1));
    end
    legend({'training data', 'promp prior', 'ground truth'});

    [visited, joint_KF] = mainloop(pmp, prm, obs);

    figure(h_error);
    plot(time(visited.k), visited.mse, strial{trial});
    % plot2svg('error')

    % return reference that is useful for ilc
    ref = prepare_reference(datac{k}.prm.ref, time, joint_KF);

    % remove path so that promp stuff does not clash with ilc
%     for k=1:numel(extra_folders)
%         rmpath([path1 '/' extra_folders{k}]);
%     end


end %

function ref = prepare_reference(ref_original, time, joint_KF)


% put reference in useful form here
% expected format of ILC usage
% 
%     param: [1x1 struct]
%      name: 'rD'
%         t: [1x151 double]
%         T: [4x4x151 double]
%         q: [151x4 double]
%       xyz: [151x3 double]
%         u: [151x4 double]
%      qdot: [151x4 double]
    
%     t = time(1:end-1);
%     tor = ref_original.t;
    
    dbg=0;
    if dbg
        figurew 'debugq1'
        plot(t,   joint_KF(1).q_mean, 'ro-');
        plot(tor, ref_original.q(:,2), 'b')

        figurew 'debugq2'
        plot(t,   joint_KF(2).q_mean, 'ro-');
        plot(tor, ref_original.q(:,3), 'b')    

        figurew 'debugq3'
        plot(t,   joint_KF(3).q_mean, 'ro-');
        plot(tor, ref_original.q(:,4), 'b')    

        figurew 'debugu1'
        plot(t,   joint_KF(4).q_mean, 'ro-');
        plot(tor, ref_original.u(:,2), 'b')    

        figurew 'debugu2'
        plot(t,   joint_KF(5).q_mean, 'ro-');
        plot(tor, ref_original.u(:,3), 'b')    

        figurew 'debugu3'
        plot(t,   joint_KF(6).q_mean, 'ro-');
        plot(tor, ref_original.u(:,4), 'b');
    end
    
    ref = ref_original;
    
    u = [];
    for j=4:6;
        u_sol = [joint_KF(j).q_mean   joint_KF(j).q_mean(end)];
        u     = [u;    interp1(time,  u_sol, ref.t)];
    end
    u = [zeros(1, numel(ref.t)); u]';
    
    if sum( isnan( u(:,end) ) ) ~= 0
        disp('nan detected when creating new ilc guess');
        disp('trying to solve it now...');
        u = [u(1:end-1,:);  zeros(1,4) ];
    end
    
    if sum( isnan( u(:,end) ) ) ~= 0
        disp('could not solve nan problem');
        disp('debug code...');
    else
        disp('NaN solved!!!!!!!!!!!');
    end
   
    ref.u_conditioned = u;
    ref.T = [];
    ref.q = [];
    ref.xyz = [];
    ref.qdot = [];
    ref.u = [];

        
end



    
%     figurew('start_ilc');
%     
% 
%     subplot(3,1,1); hold on;
%     u1_guess = interp1(sol(4).t, sol(4).q, ref.t);
%     plot(ref.t,     u1_guess, SBLUE(2));
%     plot(ref.t,   ref.u(:,2), SRED(2));
%     
%     subplot(3,1,2); hold on;
%     u2_guess = interp1(sol(5).t, sol(5).q, ref.t);
%     plot(ref.t,     u2_guess, SBLUE(2));
%     plot(ref.t,   ref.u(:,3), SRED(2));
%     
%     subplot(3,1,3); hold on;
%     u3_guess = interp1(sol(6).t, sol(6).q, ref.t);
%     plot(ref.t,     u3_guess, SBLUE(2));
%     plot(ref.t,   ref.u(:,4), SRED(2));    
% 
%     ref.u = [0.*u1_guess'  u1_guess'  u2_guess'   u3_guess'];
%     
%     
%     if sum( isnan(u1_guess) ) ~= 0
%         newIndex = ~isnan(u1_guess);
%         ref.t = ref.t(newIndex);
%         ref.u = ref.u(newIndex,:);
%         ref.q = ref.q(newIndex,:);
%         ref.qdot = ref.qdot(newIndex,:);
%         ref.xyz = ref.xyz(newIndex,:);
%     end
%     
%     























