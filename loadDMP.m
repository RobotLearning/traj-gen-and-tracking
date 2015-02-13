%% Loads DMPs from kinesthetic teach-in data

clc; clear; close all
% load the data
M = dlmread('../Desktop/okanKinestheticTeachin_20141210/unifyData65.txt');
dof = 7;
scale = 1e-3;

% extract time and the joints - last 14 variables
Mq = M(:,end-14:end);

t = Mq(:,1);
t = scale * t;
q = Mq(:,2:dof+1);
qd = Mq(:,dof+2:end);

for i = 1:dof
    joints{i} = ['joint\_', int2str(i)];
    vel{i} = ['joint\_vel\_', int2str(i)];
    %err{i} = ['err_j_', int2str(i)];
end

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

%% check whether velocity data matches to differences in one joint
% figure;
% for j = 1:dof
% subplot(7,1,j);
% plot(t,qd(:,j),t,qd2(:,j));
% legend([vel{j},'\_real'],[vel{j},'\_diff']);
% end

%% check noise spectrum 

L = 2048; %500 * ceil(length(t)/500);
NFFT = 2^nextpow2(L);
Fs = L/(t(end)-t(1));
f = Fs/2*linspace(0,1,NFFT/2+1);
tLin = linspace(t(1),t(end),L);
qLin= interp1(t,q,tLin);
qdLin = interp1(t,qd,tLin);
qSpecs = fft(qLin,NFFT)/L;
qdSpecs = fft(qdLin,NFFT)/L;

% plot single-sided amplitude spectrum.
% figure;
% for j = 1:dof
% subplot(7,2,2*j-1);
% plot(f,2*abs(qSpecs(1:NFFT/2+1,j)));
% legend(['spectra\_',joints{j}]);
% subplot(7,2,2*j);
% plot(f,2*abs(qdSpecs(1:NFFT/2+1,j)));
% legend(['spectra\_',vel{j}]);
% end

%% Filter the signals
% fft back the signals
%
% TODO: joint angles loses DC value!

cutoff = 5;
%qSpecsZeroed = [qSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
qdSpecsZeroed = [qdSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
%qFilLin = real(ifft(qSpecsZeroed,NFFT)*2*L);
qdFilLin = real(ifft(qdSpecsZeroed,NFFT)*2*L);

% interpolate back
%qFil = interp1(tLin,qFilLin,t);
qdFil = interp1(tLin,qdFilLin,t);

figure;
for j = 1:dof
%subplot(7,2,2*j-1);
%plot(t,qFil(:,j),t,q(:,j));
%legend([joints{j},'\_filt'], joints{j})
subplot(7,1,j);
plot(t,qdFil(:,j),t,qd(:,j));
legend([vel{j},'\_filt'], vel{j});
end

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
% motion takes 3 seconds so 1.5 seconds on each side
window = 3;

% take the median index of the maximum of the first 3 joints
[~,idx] = max(abs(qdFil(:,1:3)));
idx = median(idx);
idxNew = t >= t(idx)-window/2 && t <= t(idx)+window/2;
tNew = t(idxNew);
qdNew = qdFil(idxNew);

figure;
for j = 1:dof
%subplot(7,2,2*j-1);
%plot(t,qFil(:,j),t,q(:,j));
%legend([joints{j},'\_filt'], joints{j})
subplot(7,1,j);
plot(t,qdNew(:,j));
legend([vel{j},'\_cut']);
end

% feed them all to a dmp
% sample one dmp (based on a goal position) and save to txt file
% load from SL and run it
% do ILC on it

% is initial point = final point ?

