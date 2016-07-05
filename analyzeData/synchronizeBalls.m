% Synchronizes balls by trying to match their times
%
% Checks up to res accuracy between tPred and tAct 
%
function [b_sync,t_sync,idx_sync] = synchronizeBalls(tPred,tAct,ballPred,res)

    numSamp = min(length(tAct),length(tPred));
    idx = 1;
    for i = 1:numSamp
        if abs(tPred(i) - tAct(i)) < res
            b_sync(idx,:) = ballPred(i,:);
            t_sync(idx) = tPred(i);
            idx_sync(idx) = i;
            idx = idx + 1;
        end
    end

end