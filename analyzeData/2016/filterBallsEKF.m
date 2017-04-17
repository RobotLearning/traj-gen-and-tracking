% Filter balls with EKF
% Initialize using a polynomial filter

function ballEsts = filterBallsEKF(t,b,ekfFilter,spin)

% initialize the EKF filter with nonlinear ls estimate
num_init_balls = 12;
x0 = estimateInitBall(t(1:num_init_balls),b(1:num_init_balls,:)',spin);

if spin.flag
    dim = 6;
else
    dim = 3;
end

ballEsts = zeros(length(t),2*dim+1);
% x0 = estimateInitialState(ballFlightFnc,t,b,num_init_balls);
P0 = 0.01*eye(length(x0));
ekfFilter.initState(x0,P0);
ekfFilter.linearize(0.01,0);

t_last = 0.0;
for j = 1:length(t)
    diff_t = t(j) - t_last;
    %filter.linearize(dt,0);
    ekfFilter.predict(diff_t,0);
    ekfFilter.update(b(j,:)',0);
    ballEsts(j,:) = [t(j),ekfFilter.x(:)'];   
    t_last = t(j);
end

end

% % Estimate initial state using the first 12 balls and 
% % fitting a 2nd order poly to it
% function x = estimateInitialState(ballFlightFnc,t,b,numBalls)
% 
% % initialize poly filter
% dim = 3; % 2D
% order = 2; % 2nd order polynomial
% size = 12; % size of the moving average
% 
% filter = PolyFilter(dim,order,size,ballFlightFnc);
% 
% diff_t = 0;
% for i = 1:numBalls
%     filter.update(diff_t,b(i,:));
%     diff_t = t(i+1) - t(i);    
% end
% x = filter.x(:);
% 
% end