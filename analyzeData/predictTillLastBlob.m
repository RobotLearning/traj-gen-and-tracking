
function ballPred = predictTillLastBlob(filter,t1,t3,b1,b3,t_lookup,ball_lookup)

loadTennisTableValues;
% find the index at bounce
[~,idx_bounce] = min(b3(:,3));
predTime.bounce =  t3(idx_bounce) - t_lookup;
predTime.lastBlob3 = t3(end) - t_lookup;
predTime.lastBlob1 = t1(end) - t_lookup;

dt = 1/60;
initVar = 1;
filter.initState(ball_lookup(:),initVar);
filter.linearize(dt,0);
predictHorizon = predTime.lastBlob1; % only predict till last blob3
table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;
[ballPred,ballTime,numBounce,time2PassTable] = ...
            predictBallPath(dt,predictHorizon,filter,table);