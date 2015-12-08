
function [] = check_all_possible_refs(bioRob, dt, saveFigureFolderName)

    % plot all references and joint angles here
    h.r = figurew('refs');  
    view([0 -1 0]);
    drawnow; axis 'equal'; drawnow;
     set_fig_position(  [0.0448 0.413 0.272 0.541] );
    plot3(0,0,0, 'b.');
    xlim([0      0.7]);
    ylim([-0.10   0.05]);
    zlim([-0.3  0.55]);

    
    h.q =figurew('q');
    set_fig_position( [0.291 0.189 0.241 0.77])

    h.qdot =figurew('qdot');
    set_fig_position(  [0.517 0.47 0.241 0.777] );
    
    h.u =figurew('u');
    set_fig_position( [0.724 0.485 0.241 0.777] )
    

        
    
    plot_refs = {'eA', 'eB', 'eC', 'eD', 'eE'};
    style = {SRED(2), SBLUE(2), SBLACK(2), SGREEN(2), SMAG(2)};

    
    
    for r=1:length(plot_refs )

        [xz_i, Amp, tf, rev] = get_ref_arguments(plot_refs{r});
        [t, x, z] = circ(xz_i, Amp, dt, tf, rev);
        q  = do_IK(bioRob, t, x, z);
        
        u  = compute_inv_dyn(bioRob, q, dt );

        figure(h.r);
        plot3(x, 0.*z, z, style{r});
        animate_bioRob(bioRob, q, 5);
        text( xz_i(1)+0.02, 0,  xz_i(2), [ plot_refs{r} '. Amp: ' num2str(Amp(1)) ', ' num2str(Amp(2))]);

        q_total = [q(:,2:end) u(:,2:end)];
        
        figure(h.q); 
        for j=1:3
            subplot(3,1,j); hold on; 
            ylabel(['q' num2str(j) ' (rad)'])
            plot(t, q_total(:,j), style{r}); 
        end
        
        qdot_total = diff(q_total)./dt;        
        qdot_total = [qdot_total; qdot_total(end,:)];
        figure(h.qdot)
        for j=1:3
            subplot(3,1,j); hold on;
            ylabel(['qdot' num2str(j) ' (rad/s)'])
            plot(t, qdot_total (:,j), style{r}); 
        end    
        
        figure(h.u); 
        for j=1:3
            subplot(3,1,j); hold on;
            ylabel(['u' num2str(j) ' Torque'])
            plot(t, q_total(:,j+3), style{r}); 
        end

    end
    
    
    %figure(h.r); drawnow;
    %plot2svg('refs_now.svg'); cannot save 3D figure
    if ~isempty(saveFigureFolderName)
        sspng(saveFigureFolderName, []); 
    end
    

end





















