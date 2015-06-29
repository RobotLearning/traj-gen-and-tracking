%% testing striking dmp in SL


clc; clear; close all;

h = 0.002;
folder = 'home/okankoc/robolab/barrett/saveData/';
q0 = load([folder,'q0.txt']);
qf = load([folder,'qf.txt']);
W = load([folder,'w_strike.txt']);
M = dlmread([folder,'dmp_q','.txt']);
Md = dlmread([folder,'dmp_qd','.txt']);
Mdd = dlmread([folder,'dmp_qdd','.txt']);
t = h * (1:size(M,1));

% canonical system
tau = 1; %0.6 * 1/t(end);
alpha = 25;
beta = alpha/4;
ax = 1;
% number of basis functions
numbf = 50;
pat = 'd';
can = CAN(h,ax,tau,numbf,t(end),pat);

for i = 1:7
    dmp(i) = DDMP(can,alpha,beta,qf(i),[q0(i);0;0]);
    dmp(i).w = W(:,i);
end

Mmat = M;
Mdmat = Md;
Mddmat = Mdd;
figure;
for j = 1:7
    [~,qs] = dmp(j).evolve(length(t));
    Mmat(:,j) = qs(1,:)';
    Mdmat(:,j) = qs(2,:)';
    Mddmat(:,j) = qs(3,:)';
    subplot(7,3,3*j-2);
    plot(t,qs(1,:),t,M(:,j));
    legend('matlab','SL');
    subplot(7,3,3*j-1);
    plot(t,qs(2,:),t,Md(:,j));
    legend('matlab','SL');
    subplot(7,3,3*j);
    plot(t,qs(3,:),t,Mdd(:,j));
    legend('matlab','SL');
end

initializeWAM;
Q = [M,Md]';
Qmat = [Mmat,Mdmat]';
xmat = wam.kinematics(Qmat);
x = wam.kinematics(Q);

figure;
plot3(x(1,:),x(2,:),x(3,:),'--');
hold on;
plot3(xmat(1,:),xmat(2,:),xmat(3,:),'-r');
xlabel('x');ylabel('y');zlabel('z');
legend('SL','MATLAB');
axis equal;
grid on;