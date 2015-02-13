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
plot(t,q);
legend(joints);
figure;
plot(t,qd);
legend(vel);

%% check whether velocity data matches to differences in one joint
figure;
for j = 1:dof
subplot(7,1,j);
plot(t,qd(:,j),t,qd2(:,j));
legend([vel{j},'\_real'],[vel{j},'\_diff']);
end

%% check noise spectrum 
L = 500 * ceil(length(t)/500);
NFFT = 2^nextpow2(L);
Fs = L/(t(end)-t(1));
f = Fs/2*linspace(0,1,NFFT/2+1);
tLin = linspace(t(1),t(end),L);
qLin= interp1(t,q,tLin);
qdLin = interp1(t,qd,tLin);
qSpecs = fft(qLin,NFFT)/L;
qdSpecs = fft(qdLin,NFFT)/L;

% plot single-sided amplitude spectrum.
figure;
for j = 1:dof
subplot(7,2,2*j-1);
plot(f,2*abs(qSpecs(1:NFFT/2+1,j)));
legend(['spectra\_',joints{j}]);
subplot(7,2,2*j);
plot(f,2*abs(qdSpecs(1:NFFT/2+1,j)));
legend(['spectra\_',vel{j}]);
end

%% Filter the signals
% 20 Hz seems to be a good cutoff rate
% fft back the signals

qSpecsZeroed = [qSpecs(f < 20), zeros(1,sum(f >= 20))];
qdSpecsZeroed = [qdSpecs(f < 20), zeros(1,sum(f >= 20))];
qFil = ifft(qSpecsZeroed,NFFT)*L;
qdFil = ifft(qdSpecsZeroed,NFFT)*L;

figure;
for j = 1:dof
subplot(7,2,2*j-1);
plot(t,qFil(:,j));
legend([joints{j},'\_filt'])
subplot(7,2,2*j);
plot(t,qdFil(:,j));
legend([vel{j},'\_filt']);
end

% segment motion based on maximum velocity
% feed them all to a dmp
% sample one dmp (based on a goal position) and save to txt file
% load from SL and run it
% do ILC on it

% is initial point = final point ?

