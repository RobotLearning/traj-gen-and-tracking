% Filter balls with EKF
% Initialize using a polynomial filter

function ballEsts = filterBallsEKF(t,b,ekfFilter,spin)

% initialize the EKF filter with nonlinear ls estimate
num_init_balls = 5;
t(1:num_init_balls)
b(1:num_init_balls,:)'
x0 = estimateInitBall(t(1:num_init_balls),b(1:num_init_balls,:)',spin)

if spin.flag
    if spin.known
        dim = 3;
        x0 = [x0(1:3);x0(7:9)];
    else
        dim = 6;
    end
else
    dim = 3;
end

% adjust x0
% y vel seems to be estimated too high
x0(4) = 0.92*x0(4);
x0(5) = 0.88*x0(5);
x0(6) = 0.90*x0(6);

ballEsts = zeros(length(t),2*dim+1);
% x0 = estimateInitialState(ballFlightFnc,t,b,num_init_balls);
P0 = 0.01*eye(length(x0));
ekfFilter.initState(x0,P0);
ekfFilter.linearize(0.01,0);

t_last = 0.0;
for j = 1:length(t)
    diff_t = t(j) - t_last;
    %filter.linearize(dt,0);
    vel = ekfFilter.x(4:6);
    ekfFilter.predict(diff_t,0);
    % if ball bounces change spin
    if (spin.flag && spin.known && ...
            ekfFilter.x(6) > 0 &&  vel(3) < 0)
        e_t = 0.86;
        mu = 0.15;
        ballRadius = 0.02;
        xdot = reboundSpinModel([vel;spin.est],e_t,mu,ballRadius);
        spin.est = xdot(4:6);
        xLast = ekfFilter.x;
        PLast = ekfFilter.P;
        [ekfFilter,~] = initFilterEKF(spin);
        ekfFilter.initState(xLast,PLast);
        ekfFilter.linearize(0.01,0);
    end
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