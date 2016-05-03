%% Predict next ball location given current ball position and vel


function F = predictNextBall(b,C,g,data,len)

t = repmat(data(:,1)',6,1);
x = [];
for i = 1:length(len)
    vec = b((i-1)*6+1:i*6);
    x = [x,repmat(vec,1,len(i))];
end
%x = repmat(b,1,length(t));
ballPosData = data(:,2:4)';

ballFlightFnc = @(x) [x(4:6,:);ballFlightVectorModel(x(4:6,:),C,g)];
k1 = t .* ballFlightFnc(x);
x_k1 = x + k1/2;
k2 = t .* ballFlightFnc(x_k1);
x_k2 = x + k2/2;
k3 = t .* ballFlightFnc(x_k2);
x_k3 = x + k3;
k4 = t .* ballFlightFnc(x_k3);
xNext = x + (k1 + 2*k2 + 2*k3 + k4)/6;

err = xNext(1:3,:) - ballPosData;

F = diag(err'*err);

end

function xddot = ballFlightVectorModel(xdot,C,g)

v = sqrt(xdot(1,:).^2 + xdot(2,:).^2 + xdot(3,:).^2);
xddot(1,:) = -C * v .* xdot(1,:);
xddot(2,:) = -C * v .* xdot(2,:);
xddot(3,:) = g - C * v .* xdot(3,:);

end
