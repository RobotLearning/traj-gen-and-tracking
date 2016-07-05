% Calculates the prediction error (RMS) between predicted ball path
% at certain times and the actual recorded noisy trajectory [merged]
%
% Synchronizes the times first to compare the two sets
%
% NOTE: model-free nonparametric approach to get ballPred (e.g. kNN)

function err_rms = calcPredError(tPred,ballPred,tAct,ballAct)

    [ballPredSync,tPredSync,idxSync] = synchronizeBalls(tPred,tAct,ballPred,1e-2);
    differ = ballPredSync - ballAct(idxSync,:);
    err_rms = sqrt(sum(diag((differ)*(differ)'))/size(differ,1));

end