%% Loads DMPs from kinesthetic teach-in data
%

clc; clear; close all

set = 15:65; % dataset of demonstrations
dropSet = [17,20,22,23,24,29,56]; % bad recordings
for i = 1:length(dropSet)
    set = set(set ~= dropSet(i));
end

dof = 7; % seven degrees of freedom 
scale = 1e-3; % recorded in milliseconds
Q = cell(length(set),1); % concatenate to huge q and qd matrices
Qd = cell(length(set),1);
t_strike = cell(length(set),1); % different time profiles

% load the data
for i = 1:length(set)
    
    Ms = dlmread(['../Desktop/okanKinestheticTeachin_20141210/unifyData', ...
                 int2str(set(i)),'.txt']);

    % extract time and the joints - last 14 variables
    Mq = Ms(:,end-14:end);

    t = Mq(:,1);
    t = scale * t;
    q = Mq(:,2:dof+1);
    qd = Mq(:,dof+2:end);

    for j = 1:dof
        joints{j} = ['joint\_', int2str(j)];
        vel{j} = ['joint\_vel\_', int2str(j)];
        acc{j} = ['joint\_acc\_', int2str(j)];
        %err{i} = ['err_j_', int2str(i)];
    end
    
    %{
    qd2 = diff(q)./(repmat(diff(t),1,dof));
    qd2 = [qd2; qd2(end,:)];
        
    figure;
    for j = 1:dof
        subplot(7,2,2*j-1);
        plot(t,q(:,j));
        legend(joints{j});
        subplot(7,2,2*j);
        plot(t,qd(:,j));
        legend(vel{j});
    end
    %}

    %% check whether velocity data matches to differences in one joint
    %{
    figure;
    for j = 1:dof
        subplot(7,1,j);
        plot(t,qd(:,j),t,qd2(:,j));
        legend([vel{j},'\_real'],[vel{j},'\_diff']);
    end
    %}

    %% check trajectory spectrum 

    L = 2048; %500 * ceil(length(t)/500);
    NFFT = 2^nextpow2(L);
    Fs = L/(t(end)-t(1));
    f = Fs/2*linspace(0,1,NFFT/2+1);
    tLinStrike = linspace(t(1),t(end),L);
    qLin= interp1(t,q,tLinStrike);
    qdLin = interp1(t,qd,tLinStrike);
    qSpecs = fft(qLin,NFFT)/L;
    qdSpecs = fft(qdLin,NFFT)/L;

    % plot single-sided amplitude spectrum.
    %{
    figure;
    for j = 1:dof
        subplot(7,2,2*j-1);
        plot(f,2*abs(qSpecs(1:NFFT/2+1,j)));
        legend(['spectra\_',joints{j}]);
        subplot(7,2,2*j);
        plot(f,2*abs(qdSpecs(1:NFFT/2+1,j)));
        legend(['spectra\_',vel{j}]);
    end
    %}

    %% Filter the signals
    % TODO: fft back the signals problematic - some joints lose DC value!
    %

    cutoff = 50;
    %{
    qSpecsZeroed = [qSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
    qdSpecsZeroed = [qdSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
    qFilLin = real(ifft(qSpecsZeroed,NFFT)*2*L);
    qdFilLin = real(ifft(qdSpecsZeroed,NFFT)*2*L);
    %}
    w_nyquist = floor(length(t)/2);
    w_cut = cutoff/w_nyquist;
    [B,A] = butter(2,w_cut);
    qFilLin = filtfilt(B,A,qLin);
    qdFilLin = filtfilt(B,A,qdLin);
    
    % interpolate back
    qFil = interp1(tLinStrike,qFilLin,t);
    qdFil = interp1(tLinStrike,qdFilLin,t);
    
    qdd = diff(qdFil)./(repmat(diff(t),1,dof));
    qdd = [qdd; qdd(end,:)];

    %{
    figure;
    for j = 1:dof
        subplot(7,3,3*j-2);
        plot(t,qFil(:,j),t,q(:,j));
        legend([joints{j},'\_filt'], joints{j})
        subplot(7,3,3*j-1);
        plot(t,qdFil(:,j),t,qd(:,j));
        legend([vel{j},'\_filt'], vel{j});
        subplot(7,3,3*j);
        plot(t,qdd(:,j));
        legend(acc{j});
    end
    %}

    %% Segment the signals
    % segment motion based on maximum velocity
    % does not work due to artifacts!
    velMax = max(abs(qdFil));
    % cutOffVel = 0.05;
    % velThres = cutOffVel * velMax;
    % velThres = repmat(velThres,length(t),1);
    % idxMat = abs(qdFil) > velThres;
    % % find the joint which lasts the longest
    % [~,idxJoint] = max(sum(idxMat));
    % tNew = t(idxMat(:,idxJoint));
    % qdNew = qdFil(idxMat(:,idxJoint),:);


    % throw away bad examples (~5 demonstrations)
    % segment based on position and guess
    % consider the first 3 joints
    % motion takes 1.5 seconds so 0.75 seconds on each side
    windowStrike = 0.5;
    % we assume returning follows immediately after and takes about the
    % same time duration
    windowReturn = 0.5;

    % take the median index of the maximum of the first 3 joint velocities
    [~,idx] = max(abs(qdFil(:,1:3)));
    idx = median(idx);
    idxStrike = find((t >= t(idx)-windowStrike/2) & (t <= t(idx)+windowStrike/2));
    idxReturn = find((t >= t(idx)+windowReturn/2) & (t <= t(idx)+(3*windowReturn/2)));
    tStrike = t(idxStrike);
    tStrike = tStrike - tStrike(1);
    tReturn = t(idxReturn);
    tReturn = tReturn - tReturn(1);
    qddStrike = qdd(idxStrike,:);
    qddReturn = qdd(idxReturn,:);
    qdStrike = qdFil(idxStrike,:);
    qdReturn = qdFil(idxReturn,:);
    qStrike = q(idxStrike,:);
    qReturn = q(idxReturn,:);

    %{
    figure;
    for j = 1:dof
        %subplot(7,2,2*j-1);
        %plot(t,qFil(:,j),t,q(:,j));
        %legend([joints{j},'\_filt'], joints{j})
        subplot(7,2,2*j-1);
        plot(tStrike,qStrike(:,j));
        legend([joints{j},'\_strike']);
        subplot(7,2,2*j);
        plot(tStrike,qdStrike(:,j));
        legend([vel{j},'\_strike']);
    end
    figure;
    for j = 1:dof
        %subplot(7,2,2*j-1);
        %plot(t,qFil(:,j),t,q(:,j));
        %legend([joints{j},'\_filt'], joints{j})
        subplot(7,2,2*j-1);
        plot(tReturn,qReturn(:,j));
        legend([joints{j},'\_return']);
        subplot(7,2,2*j);
        plot(tReturn,qdReturn(:,j));
        legend([vel{j},'\_return']);
    end
    %}

    %close all;
    Q_strike{i} = qStrike;
    Qd_strike{i} = qdStrike;
    Qdd_strike{i} = qddStrike;
    Q_return{i} = qReturn;
    Qd_return{i} = qdReturn;
    Qdd_return{i} = qddReturn;
    t_strike{i} = tStrike;
    t_return{i} = tReturn;

end

% feed them all to a dmp
dmpStrike = trainMultiDMPs(t_strike,Q_strike,Qd_strike,Qdd_strike);
dmpReturn = trainMultiDMPs(t_return,Q_return,Qd_return,Qdd_return);

%% Write to text file

% evolve dmps and save to a text file
% put into matrix form
% t q1 qd1 q2 qd2 ... q7 qd7
Ms = [];
Mr = [];

% scale up to 500 Hz 
scale = 500/200;
tLinStrike = linspace(tStrike(1),tStrike(end),scale*length(tStrike));
tLinReturn = linspace(tReturn(1),tReturn(end),scale*length(tReturn));
dmpStrike(1).can.changeSamplingTime(0.002);
dmpReturn(1).can.changeSamplingTime(0.002);
Ms = tLinStrike(:);
Mr = tLinReturn(:);

figure;
for j = 1:dof
    [~,qs] = dmpStrike(j).evolve(length(tLinStrike));
    subplot(7,3,3*j-2);
    plot(tLinStrike,qs(1,:));
    Ms = [Ms, qs(1,:)'];
    legend([joints{j},'\_dmpStrike']);
    subplot(7,3,3*j-1);
    plot(tLinStrike,qs(2,:));
    Ms = [Ms, qs(2,:)'];
    legend([vel{j},'\_dmpStrike']);
    subplot(7,3,3*j);
    plot(tLinStrike,qs(3,:));
    legend([acc{j},'\_dmpStrike']);
end

figure;
for j = 1:dof
    [~,qr] = dmpReturn(j).evolve(length(tLinReturn));
    subplot(7,3,3*j-2);
    plot(tLinReturn,qr(1,:));
    Mr = [Mr, qr(1,:)'];
    legend([joints{j},'\_dmpReturn']);
    subplot(7,3,3*j-1);
    plot(tLinReturn,qr(2,:));
    Mr = [Mr, qr(2,:)'];
    legend([vel{j},'\_dmpReturn']);
    subplot(7,3,3*j);
    plot(tLinReturn,qr(3,:));
    legend([acc{j},'\_dmpReturn']);
end

dlmwrite('dmp_strike.txt',Ms,'delimiter','\t','precision',6);
dlmwrite('dmp_return.txt',Mr,'delimiter','\t','precision',6);
save('dmpStrike.mat','dmpStrike');
save('dmpReturn.mat','dmpReturn');

for i = 1:dof,
    Ws(:,i) = dmpStrike(i).w;
    Wr(:,i) = dmpReturn(i).w;
end
dlmwrite('w_strike.txt',Ws,'delimiter','\t','precision',6);
dlmwrite('w_return.txt',Wr,'delimiter','\t','precision',6);

%% see how well we extend to particular demonstrations

i = 10;
q_act = interp1(t_strike{i},Q_strike{i},tLinStrike,'linear','extrap');
qd_act = interp1(t_strike{i},Qd_strike{i},tLinStrike,'linear','extrap');
qdd_act = interp1(t_strike{i},Qdd_strike{i},tLinStrike,'linear','extrap');

figure;
for j = 1:dof
    dmpStrike(j).setGoal(q_act(:,j));
    dmpStrike(j).setInitState(q_act(1,j));
    dmpStrike(j).resetStates();
    [~,qs] = dmpStrike(j).evolve(length(tLinStrike));
    subplot(7,3,3*j-2);
    plot(tLinStrike,qs(1,:),tLinStrike,q_act(:,j)');
    legend([joints{j},'\_dmpStrike'],'demonstration');
    subplot(7,3,3*j-1);
    plot(tLinStrike,qs(2,:),tLinStrike,qd_act(:,j)');
    legend([vel{j},'\_dmpStrike'],'demonstration');
    subplot(7,3,3*j);
    plot(tLinStrike,qs(3,:),tLinStrike,qdd_act(:,j)');
    legend([acc{j},'\_dmpStrike'],'demonstration');
end

i = 5;
q_act = interp1(t_return{i},Q_return{i},tLinReturn,'linear','extrap');
qd_act = interp1(t_return{i},Qd_return{i},tLinReturn,'linear','extrap');
qdd_act = interp1(t_return{i},Qdd_return{i},tLinReturn,'linear','extrap');

figure;
for j = 1:dof
    dmpReturn(j).setGoal(q_act(:,j));
    dmpReturn(j).setInitState(q_act(1,j));
    dmpReturn(j).resetStates();
    [~,qs] = dmpReturn(j).evolve(length(tLinReturn));
    subplot(7,3,3*j-2);
    plot(tLinReturn,qs(1,:),tLinReturn,q_act(:,j)');
    legend([joints{j},'\_dmpReturn'],'demonstration');
    subplot(7,3,3*j-1);
    plot(tLinReturn,qs(2,:),tLinReturn,qd_act(:,j)');
    legend([vel{j},'\_dmpReturn'],'demonstration');
    subplot(7,3,3*j);
    plot(tLinReturn,qs(3,:),tLinReturn,qdd_act(:,j)');
    legend([acc{j},'\_dmpReturn'],'demonstration');
end