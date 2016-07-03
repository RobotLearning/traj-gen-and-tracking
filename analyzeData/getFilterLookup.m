
function [ball_lookup,idx_lookup,t_lookup] = getFilterLookup(ballEst)

toly = 0.4;
% load table parameters
loadTennisTableValues;

% particular rule implemented in lookup table in SL
y_center = dist_to_table - table_length/2;
idx_filter_SL = find(ballEst(:,3) > y_center & ballEst(:,3) < y_center + toly & ...
               ballEst(:,6) > 0.5,1);

ball_lookup = ballEst(idx_filter_SL,2:end);
t_lookup = ballEst(idx_filter_SL,1);
idx_lookup = idx_filter_SL;