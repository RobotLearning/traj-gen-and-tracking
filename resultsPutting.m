%% Putting results

X = 1:10;
err = zeros(10,10);

% run testPutting ten times and store the errors
for trial = 1:length(X)
    testPutting;
    err(trial,:) = ilc.error;
end

N = length(t);
% get the root mean square error from sse
RMS = sqrt(err/N);

E = std(RMS);
M = sum(RMS)/length(X);
errorbar(X,M,E/2);

xlabel('Iterations');
ylabel('RMS Error');
legend('wILC','Location','NorthEast');
title('Learning Performance in Putting');

