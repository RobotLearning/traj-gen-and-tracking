% function to find the index of the balldata closest to rebound in a robust way
function idx = findReboundIndex(ball)

    loadTennisTableValues();
    % check if there is a bounce
    tol = 5e-2;    
    idxBallBounce = [];
    y_center = dist_to_table - table_length/2;
    while isempty(idxBallBounce)
        idxBallBounce = find(ball(:,3) <= table_z + ball_radius + tol & ...
                             ball(:,3) >= table_z & ...
                             abs(ball(:,2) - y_center) < table_length/2 & ...
                             abs(ball(:,1)) < table_width/2, 1);
        tol = 1.2 * tol;
    end
    idx = idxBallBounce(1);

    
end