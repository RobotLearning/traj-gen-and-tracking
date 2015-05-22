%% Zero phase filtering

N = 101;
t = linspace(0,1,N);
signal = 5*sin(10*t) + 10*cos(5*t);
noise = 0.3 * randn(1,N);

x = signal + noise;
w = 25; % cutoff freq
w_perc = 25/floor(N/2);

[B,A] = butter(2,w_perc);
y = filtfilt(B,A,x);

% filter, reverse data, filter again, and reverse data again
y2 = filterButter2nd(x,w_perc);
y2 = y2(end:-1:1);
y2 = filterButter2nd(y2,w_perc);
y2 = y2(end:-1:1);

plot(t,x,t,y,t,y2);
legend('signal','filtfiltered signal','butterworth');