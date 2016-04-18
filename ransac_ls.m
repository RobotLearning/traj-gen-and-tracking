% RANSAC algorithm for robust estimation/fitting for linear least squares
%
% data: a 2xn dataset with n data points
% num: the number of points for fitting
% iter: the number of iterations
% threshDist: the threshold of the distances between points and the fitting line
% outlierRatio: the threshold of the number of outliers 

function betaBest = ransac_ls(x,y,num,iter,threshDist,outlierRatio)

m = size(x,2); % Dimension of fitting
n = size(x,1); % Total number of points
assert(n == size(y,1),'Number of points do not match!');
bestInlierNum = 0; % Best fitting line with largest number of inliers
betaBest = zeros(m,1);
inlierRatio = 1 - outlierRatio;

for i = 1:iter
    
    % Randomly select some points
    idx = randperm(n,num); 
    sample_x = x(idx,:);   
    sample_y = y(idx,:);
    
    beta = sample_x \ sample_y;
    distance = y - x*beta;
    
    %% Compute the inliers with distances smaller than the threshold
    inlierIdx = find(abs(distance) <= threshDist);
    inlierNum = length(inlierIdx);
    
    %% Update the number of inliers and fitting model if better model is found     
    if inlierNum >= round(inlierRatio*n) && inlierNum > bestInlierNum
        bestInlierNum = inlierNum;
        betaBest = beta;
    end
end