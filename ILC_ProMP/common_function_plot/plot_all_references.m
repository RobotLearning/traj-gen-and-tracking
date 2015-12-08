
function h = plot_all_references(dfv, bioRob, refc)


    nRef = numel(refc);
    cl  = color_finder( nRef );

    trajType = refc{1}.param.type;
    

    h.anim = figurew('bioRob'); view(dfv);
    drawnow; axis 'equal'; drawnow;
    set_fig_position(  [0.119 0.542 0.264 0.436]  );
    % Define axes limits depending on the task
    % =================================
    ylim([-0.05   0.05]);
    if strcmp(trajType, 'ellipse')
        xlim([0       0.6]);
        zlim([-0.25   0.45]);
        txti = numel(refc{1}.xyz(:,1));
    end
    if strcmp(trajType, 'p2p')
        xlim([  -.10       0.7]);
        zlim([  -0.35   0.45]);
        txti = round(numel(refc{1}.xyz(:,1))/2);
    end
    
    
    h.q = figurew('joint_position');  
    set_fig_position( [0.588 0.43 0.17 0.546]  );
    
    
    
    h.ilc = figurew('feedforward_motor_commands');
    set_fig_position([0.417 0.429 0.17 0.546]   );

    
    
    for k=1:nRef
      
        ref = refc{k};
        style = struct('Color', cl(k,:), 'LineStyle', '-', 'LineWidth', 2);
        
        names{k} = ref.name;
        
        figure(h.anim)
        plot3(ref.xyz(:,1), ref.xyz(:,2), ref.xyz(:,3), style); 
        text(ref.xyz(txti,1)+.01, ref.xyz(txti,2), ref.xyz(txti,3),...
             [ref.name])

        figure(h.q)
        for j=1:3
            subplot(3,1,j); hold on;
            plot(ref.t, r2d(ref.q(:,j+1)), style);
            xlabel('Time (s)');
            ylabel(['Joint' num2str(j+1) ' (deg)']);        
        end  
        subplot(3,1,1);
        legend(names);
        
        figure(h.ilc)        
        for j=1:3
            subplot(3,1,j); hold on;
            plot(ref.t, ref.u(:,j+1), style);
            xlabel('Time (s)');
            ylabel(['u' num2str(j+1) ]);
        end
        subplot(3,1,1);
        legend(names);        
        
    end

% Do this just to plot the robot in its initial configuration
figure(h.anim)
animate_bioRob(bioRob, refc{1}.q(1:5,:), 1);


end
    