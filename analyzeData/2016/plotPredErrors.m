% Plots the RMS prediction errors
% x axis shows the number of balls added to the filter
% y axis shows RMS error
%
% We want to see here the correction due to more balls:
% prediction error should go down as more balls are observed

function plotPredErrors(rms_errors)

    figure;
    title('Prediction error as more balls are observed');
    xlabel('Number of balls observed after initial prediction');
    ylabel('RMS of traj after observed ball');
    hold on;
    
    valid_set = [];
    for i = 1:length(rms_errors)
        if sum(isnan(rms_errors{i})) == 0
            valid_set = [valid_set,i];
        end
    end
    n_valid = length(valid_set);
        
    if iscell(rms_errors) && n_valid > 1

        % find cell with max entry
        Nmax = 1;
        for i = valid_set
            if length(rms_errors{i}) > Nmax
                Nmax = length(rms_errors{i});
            end
        end
        m = zeros(1,Nmax);
        n = zeros(1,Nmax);
        s = zeros(1,Nmax);    
        % go through each cell and collect statistics 
        % for errorplot
        for i = valid_set
            for j = 1:length(rms_errors{i})
                m(j) = m(j) + rms_errors{i}(j);
                n(j) = n(j) + 1;
                s(j) = s(j) + rms_errors{i}(j)^2;
            end
        end
        % get mean and stdev. estimates
        mean_est = m./n;
        std_est = sqrt(s./n - mean_est.^2);
        errorbar(mean_est,std_est);    
        hold off;
    else
        if iscell(rms_errors)
            rms_errors = rms_errors{1};
        end
        
        plot(rms_errors);
        hold off;
        % give info on average
        avg = sum(rms_errors)/length(rms_errors);
        fprintf('Avg RMS: %f \n', avg);
    end

end