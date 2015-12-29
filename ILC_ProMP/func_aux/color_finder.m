
function out =  color_finder( n )

    out = [];
    if exist('n', 'var')
        %out = distinguishable_colors( n );
        func = @(x) colorspace('RGB->Lab',x);
        out = distinguishable_colors(n,'w',func);
    else
        standard_display;
    end
    
    
end

function [] = standard_display()

    rgb.r = [0:0.25:1];
    rgb.g = [0:0.25:1];
    rgb.b = [0:0.25:1];

    figurew('rgb_pallete'); hold on; grid on;
    counter= 0;

    for ir = 1:length(rgb.r)
        for ig = 1:length(rgb.g)
            for ib = 1:length(rgb.b)
               rgb.rgb = [rgb.r(ir)  rgb.g(ig) rgb.b(ib)];
               fprintf('Counter %g, rgb [%g %g %g]\n', counter, rgb.rgb);
               style = struct('Color', rgb.rgb, 'LineStyle', 'none', ...
                              'Marker', 'o', 'MarkerFaceColor', ...
                               rgb.rgb, 'LineWidth', 2, 'MarkerSize', 40);
               plot(counter, 0, style) 
               counter = counter+1;
            end    
        end    
    end

    xlim([1 counter])
    ylim([-1  1])
    set_fig_position( [0.0682 0.611 0.882 0.281]   );

end