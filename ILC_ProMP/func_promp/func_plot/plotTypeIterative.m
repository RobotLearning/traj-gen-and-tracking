
function [h] = plotTypeIterative(h, time, nTraj, data, prm)


    nSize = nTraj;
    nJntA = prm.nJointA;
    nJntB = prm.nJointB;

    
    %time = [1:nSize-1];
    
    color_ = 'r';
    
    
    time = time(1:end-1);
    % ===================================
    % Plot position of joints
    % ===================================
    figure(h.qA);
    
    for j=1:nJntA
        subplot(1,nJntA,j); 
        Hq{j}      = shadedErrorBar(time, data(j).q_mean , 2*sqrt( data(j).q_Sigma_ii ) , ...
                                  {color_, 'LineWidth', 1'}, 1);
        delete(Hq{j}.edge);
        h.qshade = [h.qshade    Hq{j}.patch   Hq{j}.mainLine];
    end
    legend({'training data', 'promp prior', 'desired trask', 'observations', 'promp posterior'});

    
    % ===================================
    % Plot angular velocity of joints
    % ===================================
    figure(h.qAdot);
    
    for j=1:nJntA
        subplot(1,nJntA,j);
        Hq{j}      = shadedErrorBar(time, data(j).qdot_mean , 2*sqrt( data(j).qdot_Sigma_ii ) , ...
                                  {color_, 'LineWidth', 1'}, 1);
        delete(Hq{j}.edge);
        h.qshade = [h.qshade    Hq{j}.patch   Hq{j}.mainLine];
    end
    legend({'training data', 'promp prior', 'desired task', 'observations', 'promp posterior'});

    
    % ===================================
    % Plot estimated commands
    % ===================================
    figure(h.qB);
    for j=1:nJntB
        subplot(1,nJntB,j);
        Hq{j}      = shadedErrorBar(time, ...
                               data(j+nJntA).q_mean , ...
                               2*sqrt( data(j+nJntA).q_Sigma_ii ) , ...
                                  {color_, 'LineWidth', 1'}, 1);
        delete(Hq{j}.edge);
        h.qshade = [h.qshade    Hq{j}.patch   Hq{j}.mainLine];
    end
    legend({'training data', 'promp prior', 'ground truth', 'promp posterior (ILC guess)'});

   
    
end
