% Generate fitted ball path

function xNext = fitBallPath(b0,C,g,t)


t = repmat(t',6,1);
x = repmat(b0,1,length(t));

ballFlightFnc = @(x) [x(4:6,:);ballFlightVectorModel(x(4:6,:),C,g)];
k1 = t .* ballFlightFnc(x);
x_k1 = x + k1/2;
k2 = t .* ballFlightFnc(x_k1);
x_k2 = x + k2/2;
k3 = t .* ballFlightFnc(x_k2);
x_k3 = x + k3;
k4 = t .* ballFlightFnc(x_k3);
xNext = x + (k1 + 2*k2 + 2*k3 + k4)/6;

end

function xddot = ballFlightVectorModel(xdot,C,g)

v = sqrt(xdot(1,:).^2 + xdot(2,:).^2 + xdot(3,:).^2);
xddot(1,:) = -C * v .* xdot(1,:);
xddot(2,:) = -C * v .* xdot(2,:);
xddot(3,:) = g - C * v .* xdot(3,:);

end