function out = style(varargin)
% function out = style(varargin)
%   out = style(color, marker, LineWidth, LineStyle, sizeBall)
%   ex1:  style('r',   'o')
%   ex2:  style('r',   [], [], '--', 10)
%
% Default values
%     color       blue
%     marker      'none'
%     lineweight  1.5
%     linestyle   '-'
%     sizeBall    7

%           out = style()
%           out = style('r', '--')
%           out = style('r', '--', 3)

    default{1} = 'b';
    default{2} = 'none';
    default{3} = 1.5;
    default{4} = '-';
    default{5} = 7;
    
    for k=1:numel(varargin)
        try varargin{k};
            if ~isempty(varargin{k})
                default{k} = varargin{k};
            end
        end
    end
     
    out = struct('Color', default{1}, ...
                 'Marker' , default{2},...
                 'LineWidth', default{3},...
                 'LineStyle', default{4}, ...
                 'MarkerSize', default{5},...
                 'MarkerFaceColor', default{1} );
               
end






















