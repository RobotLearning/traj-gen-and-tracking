% This function prunes the really bad ball data
% Removes the unreliable - marked as 0 in the text file - ball data
%
% b1 is the pruned data from cameras 1 and 2
% b3 is the pruned data from cameras 3 and 4

function [b1,b3,filterSL] = pruneBallData(t,B)

% load table parameters
loadTennisTableValues;

% this is to throw away really bad observations as to not confuse the
% filter
zMax = 0.5;

% camera 1 time indices for reliable observations
idx1 = find(B(:,2) == 1 & B(:,5) >= table_z & B(:,5) < zMax);
t1 = t(idx1);
% order the ball data w.r.t. time and throw away same values
b1 = [t1,B(idx1,3:5)];
b1 = sortrows(b1);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos1 = []; % index for diff ball positions
for i = 1:size(b1,1)-1
    if norm(b1(i,2:4) - b1(i+1,2:4)) > tol
        idxDiffBallPos1(j) = i;
        j = j+1;
    end
end
b1 = b1(idxDiffBallPos1,:);

% camera 3 time indices for reliable observations
idx3 = find(B(:,7) == 1 & B(:,10) >= table_z & B(:,10) < zMax);
t3 = t(idx3);
% order the ball data w.r.t. time and throw away same values
b3 = [t3,B(idx3,8:10)];
b3 = sortrows(b3);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos3 = []; % index for diff ball positions
for i = 1:size(b3,1)-1
    if norm(b3(i,2:4) - b3(i+1,2:4)) > tol
        idxDiffBallPos3(j) = i;
        j = j+1;
    end
end
b3 = b3(idxDiffBallPos3,:);

% if available get the filter data
try
    b3est = [t3,B(idx3,12:17)];
    b3est = sortrows(b3est);
    filterSL = b3est(idxDiffBallPos3,:);
catch 
    warning('SL filter is not available!');
    filterSL = [];
end