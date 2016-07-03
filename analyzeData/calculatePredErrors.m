% Calculate the prediction errors using ball model (method of filter class)
% BallEsts are provided by (possibly another) filter estimates ballEsts
% 
% We test here the reduction of error as we get more ball data
% we test with the ball estimates (ballEsts) indexed from idx_start
% to idx_end 

function RMS_pred = calculatePredErrors(filter,ballEsts,idx_start,idx_end,t,b)

N_updates = idx_end - idx_start;
RMS_pred = zeros(N_updates,1);

for i = 1:N_updates
    t_now = ballEsts(idx_start+i,1);
    ball_now = ballEsts(idx_start+i,2:end);
    ballPredNew = predictTillLastBlob(filter,t,t_now,ball_now);
    x_pred = ballPredNew(1:3,:)';
    x_meas = b(idx_start+i+1:end,:);
    differ = x_pred - x_meas;
    N = size(x_pred,1);
    RMS_pred(i) = sqrt(sum(diag((differ)*(differ)'))/N);
end

disp('Plotting last prediction result also');
plotPredictionResult([],b,ballPredNew);