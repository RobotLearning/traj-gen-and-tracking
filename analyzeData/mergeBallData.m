% Merge ball data from cameras 1 and 3
% Assumes that the data is preprocessed, outliers are removed etc.

function [tMerge,ballMerge] = mergeBallData(t1,b1,t3,b3)

% find the index at which t3 stops and t1 continues
idx1butNot3 = t1 > t3(end);
ballMerge = [b3; b1(idx1butNot3,:)];
tMerge = [t3; t1(idx1butNot3)];

end