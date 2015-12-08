function [] = set_fig_position( b )
% b = [1x4]
% is a vector that is obtained by using the command getBoxSize for the
% current active figure.

    set(gcf,'Units','normalized');
    set(gcf, 'OuterPosition' , b );

end