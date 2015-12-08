function [style] = SGRAY(lw, strength)



    if ~exist('lw')
        lw = 1;
    end
    
    if exist('strength')
        strength = [strength strength strength];
    else
        strength = [0.7 0.7 0.7];
    end

    style = struct('Color', strength, 'LineStyle', '-', 'LineWidth', lw);
    
end