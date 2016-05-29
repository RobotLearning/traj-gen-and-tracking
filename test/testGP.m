%% Test Gaussian Process (GP) regression

clc; clear; close all;

%% Simple regression
hp.type = 'squared exponential ard';
hp.l = 0.25;
hp.scale = 1;
hp.noise.var = 0.01;

n = 20;
mu = 0;
x = randn(n,1);
l = hp.l;
s = hp.scale;
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = zeros(n,n);
for i = 1:n
    for j = 1:n
        Sigma(i,j) = kern(x(i),x(j));
    end
end
y = mu + chol(Sigma + hp.noise.var*eye(n)) * x;
gp = GP(hp,x',y);

% figure(1)
% set(gca, 'FontSize', 24)
% plot(x, y, '+', 'MarkerSize', 12)
% %axis([-1.9 1.9 -0.9 3.9])
% grid on
% xlabel('input, x')
% ylabel('output, y')

z = linspace(-2, 2, 100);
for i = 1:length(z)
    [mean, var] = gp.predict(z(i));
    m(i) = mean;
    s2(i) = var;
end

figure(1)
title('GP Regression');
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([z; flip(z,1)], f, [7 7 7]/8);
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
%axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
hold off;

%% Hyperparameter estimation

gp.fitHP(hp);
z = linspace(-2, 2, 100);
for i = 1:length(z)
    [mean2, var2] = gp.predict(z(i));
    m(i) = mean2;
    s2(i) = var2;
end

figure(2)
title('GP regression + HP estimation');
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([z; flip(z,1)], f, [7 7 7]/8);
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
%axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
hold off;