
function RMSE_ = RMSE(tExp, xExp, tRef, xRef)
% RMSE_ = RMSE(tExp, xExp, tRef, xRef)
% 
% This is a new code that calculates RMSE with an arbitrary J number of inputs.
%
% INPUT
%   tExp:   [1xN]: raw data from experiment (full range)
%   xExp:   [JxN]: raw data from experiment (full range)
%   tRef:   [1xM]: reference time (be sure range is less than tExp so that interpolation 
%                  will never return NaN)
%   xRef:   [JxM]: reference position
% 
% OUTPUT
%   RMSE in each of the inputs.
%
% *obs: for getting the RMSE as a distance between the points on the first and second line
%       sqrt( RMSE(1)^2 + RMSE(2)^2 )
%       If RMSE(1) is the error in x direciton
%          RMSE(2) is the error in y direciton
%          sqrt( RMSE(1)^2 + RMSE(2)^2 ) is the distance error 
%

% check conditions

    [lines1, col1] = size(xExp);
    [lines2, col2] = size(xRef);
    
    if lines1 ~= lines2
        error('xExp and xRef have different number of inputs.');
    end
    

    % put experimental values in the same clock as the tRef
    xquery = [];
    for j = 1:lines1
        xquery = [xquery; interp1(tExp, xExp(j,:), tRef)];
    end
    
    % cut the parts that could not be interpolated
    xquery = xquery(:,   ~isnan(xquery(1,:))   );
    
    ir = 1:numel(  xquery(1,:)  );
    
    % tracking error
    xError = [];
    for j = 1:lines1
        xError = [xError; xRef(j,ir) - xquery(j,ir)];
    end
    
    N = ir(end);
    
    % MSE in each direction
    RMSE_ = [];
    for j = 1:lines1
        RMSE_ = [RMSE_; sqrt( sum(xError(j,:).^2 )/ N)];
    end
    
    

end