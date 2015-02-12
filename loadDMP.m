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

%idxStart = 900;
%idxEnd = 1200;
idxStart = 1;
idxEnd = length(t);
tCut = t(idxStart:idxEnd);
qCut = q(idxStart:idxEnd,:);
qdCut = qd(idxStart:idxEnd,:);

qdCut2 = diff(qCut)./(repmat(diff(tCut),1,dof));
qdCut2 = [qdCut2; qdCut2(end,:)];

figure;
plot(tCut,qCut);
legend(joints);
figure;
plot(tCut,qdCut);
legend(vel);

% check whether velocity data matches to differences in one joint
figure;
j = 5;
plot(tCut,qdCut(:,j),tCut,qdCut2(:,j));
legend([vel{1},'\_real'],[vel{1},'\_diff']);

% check noise spectrum
L = 500 * ceil((idxEnd-idxStart)/500);
Fs = L/(tCut(end)-tCut(1));
tLin = linspace(tCut(1),tCut(end),L);
qdLin = interp1(tCut,qdCut(:,j),tLin);
NFFT = 2^nextpow2(L);
qdSpecs = fft(qdLin,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% plot single-sided amplitude spectrum.
figure;
plot(f,2*abs(qdSpecs(1:NFFT/2+1)))
title('Single-Sided Amplitude Spectrum of qd(t)')
xlabel('Frequency (Hz)')
ylabel('|qd(f)|')

% segment motion based on maximum velocity
% cut off the rest and filter
% 40 Hz seems to be a good cutoff rate
% feed them all to a dmp
% sample one dmp (based on a goal position) and save to txt file
% load from SL and run it
% do ILC on it

% is initial point = final point ?

