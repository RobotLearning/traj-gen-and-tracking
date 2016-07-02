% Merge ball data from cameras 1 and 3
% Assumes that the data is preprocessed, outliers are removed etc.

function ballMerge = mergeBallData(t1,b1,t3,b3)

% find the index at which t3 stops and t1 continues
idx1butNot3 = length(t3) + 1: length(t1);
ballMerge = [b3; b1(idx1butNot3,:)];

end