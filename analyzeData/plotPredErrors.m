% Plots the RMS prediction errors
% x axis shows the number of balls added to the filter
% y axis shows RMS error
%
% We want to see here the correction due to more balls:
% prediction error should go down as more balls are observed

function plotPredErrors(rms_errors)

figure;
plot(rms_errors);
title('Prediction error as more balls are observed');
xlabel('Number of balls observed after initial prediction');
ylabel('RMS of traj after observed ball');

% give info on average
avg = sum(rms_errors)/length(rms_errors);
fprintf('Avg RMS: %f \n', avg);


end