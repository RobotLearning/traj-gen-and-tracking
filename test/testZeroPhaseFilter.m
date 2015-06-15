%% Zero phase filtering

N = 101;
t = linspace(0,1,N);
signal = 5*sin(10*t) + 10*cos(5*t);
noise = 0.3 * randn(1,N);

x = signal + noise;
w = 25; % cutoff freq
w_perc = 25/floor(N/2);

% filter, reverse data, filter again, and reverse data again
y1 = filterButter2nd(x,w_perc);
y1 = y1(end:-1:1);
y1 = filterButter2nd(y1,w_perc);
y1 = y1(end:-1:1);

if exist('filtfilt')
    [B,A] = butter(2,w_perc);
    y2 = filtfilt(B,A,x);
    plot(t,x,t,y1,t,y2);
    legend('signal','butterworth twice','filtfilt');
else
    plot(t,x,'-',t,y1,'--');
    legend('signal','butterworth twice');
end