function [style] = SBLUEBALL(varargin)
    if isempty(varargin)
        sizeBall = 8;
    else
        sizeBall = varargin{1};
    end
    style = struct('Color', 'b', 'LineStyle', 'none', 'LineWidth', 2, 'Marker' , 'o',...
                   'MarkerFaceColor','b', 'MarkerSize', sizeBall);
end