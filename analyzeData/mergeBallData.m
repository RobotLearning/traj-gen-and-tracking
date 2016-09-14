% Merge ball data from cameras 1 and 3
% Subtracts the constant offset from camera 1 before merging
% Assumes that the data is preprocessed, outliers are removed etc.

function [tMerge,ballMerge,b1,offset] = mergeBallData(t1,b1,t3,b3)

% find the index at which t3 stops and t1 continues
idx1butNot3 = t1 > t3(end);
% subtract a constant offset from camera 1 data
offset = calculateOffset(t1,b1,t3,b3);
N = size(b1,1);
b1 = b1 - repmat(offset,N,1);
ballMerge = [b3; b1(idx1butNot3,:)];
tMerge = [t3; t1(idx1butNot3)];
% time starts from 0
tMerge = tMerge - tMerge(1);

end