% Polynomial filter as implemented in SL
% estimates 2nd order polynomial with moving average over 12 balls 

function ballEsts = filterBallsPoly(ballFlightFnc,t,b)

% initialize poly filter
dim = 3; % 2D
order = 2; % 2nd order polynomial
size = 12; % size of the moving average

filter = PolyFilter(dim,order,size,ballFlightFnc);

diff_t = 0;
ballEsts = zeros(length(t),2*dim+1);
for i = 1:length(t)-1
    filter.update(diff_t,b(i,:));
    diff_t = t(i+1) - t(i);    
    ballEsts(i,:) = [t(i),filter.x(:)'];
end
filter.update(diff_t,b(end,:));
ballEsts(end,:) = [t(end),filter.x(:)'];