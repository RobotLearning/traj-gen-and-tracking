%% Predict next ball location given current ball position and vel


function F = predictNextBall(b,C,g,data,len,spin)

if spin.flag
    dim = 12;
    Clift = spin.Clift;
    ballFlightFnc = @(x) [x(7:12,:);ballSpinModel(x(7:12,:),C,Clift,g)];
else
    dim = 6;
    ballFlightFnc = @(x) [x(4:6,:);ballFlightModel(x(4:6,:),C,g)];
end

t = repmat(data(:,1)',dim,1);
x = [];
for i = 1:length(len)
    vec = b((i-1)*dim+1:i*dim);
    x = [x,repmat(vec,1,len(i))];
end
%x = repmat(b,1,length(t));
ballPosData = data(:,2:4)';

if len == 1
    x = x';
end

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
