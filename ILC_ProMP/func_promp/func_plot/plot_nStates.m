function [ h ] = plot_nStates(time, datafull, prm, optionalName)


    data = datafull.work_format;
    
    nSize = length(data(1).q_mean);
    Tf = time(end);
    nJntA = prm.nJointA;
    nJntB = prm.nJointB; 
    
    facx=0.75;
    
%% create figure

    if exist('optionalName', 'var')
        handle_qA = figurew(['q_' datafull.taskA(1:end-1) '_' optionalName]);    
    else
        handle_qA = figurew(['q_' datafull.taskA(1:end-1)]);
    end
    
    set_fig_position( [0.304 0.667 facx*0.634 0.308]   ); 
    for j=1:nJntA
        subplot(1,nJntA,j); grid on; hold on;
        ylabel(['q' num2str(j) ' (rad)']);
        xlabel 'Time (s)';
        xlim([0 Tf]);
    end
    
    if exist('optionalName', 'var')
        handle_qAdot = figurew(['qdot_' datafull.taskA(1:end-1) '_' optionalName]);    
    else
        handle_qAdot = figurew(['qdot_' datafull.taskA(1:end-1)]);
    end    

    set_fig_position( [0.304 0.667-0.308  facx*0.634 0.308]   ); 
    for j=1:nJntA
        subplot(1,nJntA,j); grid on; hold on;
        ylabel(['qdot' num2str(j) ' (rad/s)']);
        xlabel 'Time (s)';
        xlim([0 Tf]);
    end
    
    if exist('optionalName', 'var')
        handle_qB = figurew(['u_' datafull.taskA(1:end-1) '_' optionalName]);    
    else
        handle_qB = figurew(['u_' datafull.taskA(1:end-1)]);
    end        

    set_fig_position( [0.304 0.667-2*0.308 facx*0.634 0.308]  );
    for j=1:nJntB
        subplot(1,nJntB,j); grid on; hold on;
        ylabel(['u' num2str(j) ' (Torque)']);
        xlabel 'Time (s)';
        xlim([0 Tf]);
    end
    
    % handle for figures
    h.qA        = handle_qA;
    h.qAdot     = handle_qAdot;    
    h.qB        = handle_qB;
    h.qshade    = []; %Hq{j};
    h.qdot      = [];
    h.qdotshade = []; %Hqdot{j};  
    h.obs =[];
    
    
    %% plot all demonstrations as thing grey lines
    if prm.plot.trainData
        figure(h.qA);
        for j=1:nJntA
            subplot(1, nJntA, j); 
            plot(time, data(j).q(1:end-1,:)', SGRAY(.15, .8), 'HandleVisibility','off');
            plot(time, data(j).q(end,:)', SGRAY(.15, .8));
        end
        legend({'training data'});
        
        figure(h.qAdot);
        for j=1:nJntA
            subplot(1, nJntA, j); 
            plot(time, data(j).qdot(1:end-1,:), SGRAY(.15, .8), 'HandleVisibility','off');
            plot(time, data(j).qdot(end,:)', SGRAY(.15, .8));
        end
        legend({'training data'});
        
        figure(h.qB);
        for j=1:nJntB
            subplot(1, nJntB, j); 
            plot(time, data(j+nJntA).q(1:end-1,:)', SGRAY(.15, .8), 'HandleVisibility','off');
            plot(time, data(j+nJntA).q(end,:)', SGRAY(.15, .8));
        end
        legend({'training data'});
    end

    %%
    % ===================================
    % Plot position of joints
    % ===================================
    figure(h.qA);
    for j=1:nJntA
        subplot(1,nJntA,j); 
        Hq{j}      = shadedErrorBar(time, data(j).q_mean,...
                                    2*sqrt( data(j).q_var ) , ...
                                    {'g', 'LineWidth', 2'}, 1);
        delete(Hq{j}.edge);
        delete(Hq{j}.mainLine);
    end
    legend({'training data', 'promp prior'});    
    
    % ===================================
    % Plot angular velocity of joints
    % ===================================
    figure(h.qAdot);
    for j=1:nJntA
        subplot(1,nJntA,j); 
        Hq{j}      = shadedErrorBar(time, data(j).qdot_mean , 2*sqrt( data(j).qdot_var ) , ...
                                  {'g', 'LineWidth', 2'}, 1);
        delete(Hq{j}.edge);
        delete(Hq{j}.mainLine);
    end
    legend({'training data', 'promp prior'});
        
    figure(h.qB);
    for j=1:nJntB
        subplot(1,nJntB,j); 
        Hq{j}      = shadedErrorBar(time, ...
                                    data(j+nJntA).q_mean , ...
                                    2*sqrt( data(j+nJntA).q_var ) , ...
                                    {'g', 'LineWidth', 2'}, 1);
        delete(Hq{j}.edge);
        delete(Hq{j}.mainLine);
    end
    legend({'training data', 'promp prior'});
    
    if prm.plot.sameYscale % scale y axis 
        figure(h.qA);
        for j=1:nJntA
            subplot(1,nJntA,j);
            ylim(prm.range(j,:));
        end
        figure(h.qB);
        for j=1:nJntB
            subplot(1,nJntB,j);
            ylim(prm.range(j+nJntA,:));
        end
    end
    



end



