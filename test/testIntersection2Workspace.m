% Test intersection with workspace

% create 2d points
n = 10;
X = 10 * rand(n,2);

D = delaunay(X(:,1),X(:,2));
trimesh(D,X(:,1),X(:,2))

C = convhull(X(:,1),X(:,2));
P = X(C(1:end-1),:);

% find if there is an intersection

t = 0:0.25:10;
x = randn*t + rand;

intersect = false;
p = [t(:),x(:)];
for i = 1:length(p)
    pt = p(i,:);
    intersect = checkIfPointIsInsidePolygon(P,pt);
    if intersect
        break;
    end
end

intersect

hold on;
plot(t,x)

hold off;