%% Tests kernel density estimation using a Gaussian mixture

% create a 2D gaussian mixture
N = 500;
d = 2;
%mu1 = [0,0]; 
%mu2 = [3.5 0];
mixture = @(x) (exp(-0.5 * (x(1)^2 + x(2)^2)) + ...
               exp(-0.5 * ((x(1)-3.5)^2 + x(2)^2))) / (4*pi);

% generate synthetic data
data = [randn(500,2);
        randn(500,1) + 3.5, randn(500,1)];

density = kernelDensityEstimation(data');

gridline = -2:0.25:6;
[X,Y] = meshgrid(gridline);
Z = zeros(size(X));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        pt = [X(i,j);Y(i,j)];
        Z(i,j) = density(pt);
    end
end

% plot the data and the density estimate
contour3(X,Y,Z,50), hold on
plot(data(:,1),data(:,2),'r.','MarkerSize',5)

% integrate density
area = 0.25^2 * sum(sum(Z))