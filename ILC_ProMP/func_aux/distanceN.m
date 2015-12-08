function dist = distanceN(queryPoint, path)
% dist = distanceN(queryPoint, path)
%   Compute the euclidean distance between a nStates dimensional query
%   point and all the points in a path.
% 
%   INPUT
%     queryPoint = [1 x nStates]. Ex = [2 3 6]
%     path = [nSteps x nStates] . Ex = [1 2 3; 4 5 6; 7 8 9; 1 2 3];
%
%   OUTPUT
%     dist = [nSteps x 1];

    dist = sqrt( sum(  bsxfun(@minus, queryPoint, path).^2, 2 ) );



end

%t = bsxfun(@minus, queryPoint, path)