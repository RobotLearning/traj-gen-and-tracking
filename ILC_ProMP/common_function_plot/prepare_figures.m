
function h = prepare_figures(dfv, ref, user_color, h)

    if isempty(user_color)
        style = struct('Color', 'b', 'LineStyle', '-', 'LineWidth', 2);
    else
        style = struct('Color', user_color, 'LineStyle', '-', 'LineWidth', 2);
    end

    flag_create_fig=0;
    if isempty(h) % only create these figures once.
        flag_create_fig = 1;
    end
        
    try close(h.time); end
    h.time = figurew('joint_position');  
    set_fig_position( [0.588 0.43 0.17 0.546]  );
    for j=1:3
        subplot(3,1,j); hold on;
        plot(ref.t, ref.q(:,j+1), style);
        xlabel('Time (s)');
        ylabel(['Joint ' num2str(j+1) '(rad)']);        
    end
   
    try close(h.cart); end
    h.cart = figurew('cartesian_position');  
    set_fig_position( [0.512 0.0315 0.219 0.406]   ); 
    drawnow; axis 'equal'; drawnow;
    plot(ref.xyz(:,1), ref.xyz(:,3), style);
    xlabel('x (m)');     ylabel('y (m)');

    % =================================
    % Make animation plot
    % =================================
    try close(h.anim); end
    h.anim = figurew('bioRob'); view(dfv);
    drawnow; axis 'equal'; drawnow;
    set_fig_position(  [0.119 0.542 0.264 0.436]  );
    plot3(ref.xyz(:,1), ref.xyz(:,2), ref.xyz(:,3), style);
    
    % Define axes limits depending on the task
    % =================================
    ylim([-0.05   0.05]);
    if strcmp(ref.param.type, 'ellipse')
        xlim([0       0.7]);
        zlim([-0.10   0.55]);
    end
    if strcmp(ref.param.type, 'p2p')
        xlim([  -.10       0.7]);
        zlim([-0.25   0.55]);
    end
    
    try close(h.ilc); end
    h.ilc = figurew('feedforward_motor_commands');
    set_fig_position([0.417 0.429 0.17 0.546]   );
    for j=1:3
        subplot(3,1,j); hold on;
        plot(ref.t, ref.u(:,j+1), style);
        xlabel('Time (s)');
        ylabel(['u ' num2str(j+1) ]);
    end
    
    try close(h.error); end
    h.error = figurew('joint_error');
    set_fig_position([0.76 0.431 0.17 0.546] );
    for j=1:3
        subplot(3,1,j); hold on;
        xlabel('Time (s)');
        ylabel(['Error ' num2str(j+1) ' (rad)']);
    end
        
    if flag_create_fig
        h.cart_rmse = figurew('cart_rmse');
        set_fig_position( [0.731 0.0111 0.217 0.345] );
        xlabel('Iteration number');
        ylabel(['RMS XY error (m)']);
    end
    
end


    