% plot results iteratively on top of nominal trajectory and inputs
function ilc_results(t,x_iter,u_app,iter)

figure(1);
plot(x_iter(1,:),x_iter(2,:),'-b');
hold on;

pause(0.5);

figure(2);
subplot(2,1,1)
plot(t(1:end-1),u_app(1,:,iter), '-b'); 
hold on;

subplot(2,1,2)
plot(t(1:end-1),u_app(2,:,iter), '-b');
hold on;