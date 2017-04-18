%% Plot trial data (blob data from cameras 1-4)

function plot_trial(b1,b3,ball_filt)

loadTennisTableValues;

figure;
try
    s1 = scatter3(b1(:,1),b1(:,2),b1(:,3),'r');
    s1.MarkerEdgeColor = s1.CData; % due to a bug in MATLAB R2015b
catch
    warning('No ball data provided from camera 1');
    scatter3(0.0,-2.5,-0.0,'r');
end
hold on;
s3 = scatter3(b3(:,1),b3(:,2),b3(:,3),'b');
%predColor = [0.200 0.200 0.200];
try 
    sP = scatter3(ball_filt(:,1),ball_filt(:,2),ball_filt(:,3),'k');
    sP.MarkerEdgeColor = sP.CData;
    legend('cam1','cam3','filter');
catch
    warning('Filtered values not supplied!');
    legend('cam1','cam3');
end

title('Ball trajectory in trial data');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);

s3.MarkerEdgeColor = s3.CData;
%sP.MarkerEdgeColor = sP.CData;
legend('cam1','cam3'); %'filter'
hold off;