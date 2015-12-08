
function [h] = plot_observations(h, time, nJointA, obs, visited)


    styleObs  = struct('Color', 'b', 'LineStyle', 'none', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','auto', 'MarkerSize', 10);
               
               
    try delete(h.obs); end
    
    % plot the observation at the current step
    if ~isempty(obs.joint)
        figure( h.qA );
        for jj = 1:nJointA
            subplot(1, nJointA ,jj);
            h.obs =[h.obs ...
                         plot( time(visited.k), visited.q(jj,:), styleObs)]; 
        end
    end
    legend({'training data', 'promp prior', 'ground truth',  'observations'});

    if ~isempty(obs.jointvel)
        figure( h.qAdot );
        for jj = 1:nJointA
            subplot(1, nJointA ,jj);
            h.obs = [h.obs ...
                         plot( time(visited.k), visited.qdot(jj,:), styleObs)]; 
        end
    end
    legend({'training data', 'promp prior', 'ground truth',  'observations'});

end