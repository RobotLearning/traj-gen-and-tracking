%% Putting Simulation results to put to paper 

clc; clear; close all;
% trying different algorithms ten time

algs = {'bILC','aILC','mILC'};
numAlgs = length(algs);
numTrials = 10;
numEpisodes = 10;
trials = 1:numTrials;
err = zeros(numAlgs,numTrials,numEpisodes);

% run testPutting ten times and store the errors
for alg = 1:numAlgs
    for trial = 1:numTrials
        ilcClass = algs{alg};
        testPutting;
        err(alg,trial,:) = ilc.error;
    end
end

N = length(t);
% get the root mean square error from sse
%RMS = sqrt(err/N);

figure(1);
for i = 1:numAlgs
    RMS = squeeze(err(i,:,:));
    E = std(RMS);
    M = sum(RMS)/length(trials);
    errorbar(trials,M,min(M,E),E);
    hold on;
end

xlabel('Iterations');
ylabel('RMS Error');
legend('bILC','aILC','mILC');
title('Learning Performance in Putting');

