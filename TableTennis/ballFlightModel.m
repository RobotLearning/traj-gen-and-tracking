%% Nonlinear ballistic flight model involving air drag
% no spin

function xddot = ballFlightModel(xdot,C,g)

v = sqrt(xdot(1)^2 + xdot(2)^2 + xdot(3)^2);
xddot(1) = -C * v * xdot(1);
xddot(2) = -C * v * xdot(2);
xddot(3) = g - C * v * xdot(3);

xddot = xddot(:);
end