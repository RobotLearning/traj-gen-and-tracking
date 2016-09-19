% Returns the SL filter estimate at lookup time
% 
function [ball_lookup,idx_lookup,t_lookup] = getFilterLookup(ballEst,spin)

% load table parameters
loadTennisTableValues;

% particular rule implemented in lookup table in SL
y_center = dist_to_table - table_length/2;

if spin.flag && ~spin.known
    
    idx_filter = find(ballEst(:,3) > y_center & ballEst(:,9) > 0.5,1);
    ball_lookup = ballEst(idx_filter,2:end);
    t_lookup = ballEst(idx_filter,1);
    idx_lookup = idx_filter;
    
else
    idx_filter_SL = find(ballEst(:,3) > y_center & ballEst(:,6) > 0.5,1);
    ball_lookup = ballEst(idx_filter_SL,2:end);
    t_lookup = ballEst(idx_filter_SL,1);
    idx_lookup = idx_filter_SL;
end