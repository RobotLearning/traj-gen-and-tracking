% Basic Iterative Learning Control using PD-type input update

classdef bILC < ILC
    
    % fields common to all ILCs (defined in abstract ILC class)
    properties
         
        % number of total episodes so far
        episode
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        error
        
        % ILC's Last input sequence
        u_last
    end
    
    methods
        
        function obj = bILC(trj)
                        
            obj.episode = 0;
            obj.color = 'b';
            obj.name = 'Arimoto type ILC';
            obj.error = 0;
            
            N = trj.N - 1;
            obj.u_last = trj.unom(:,1:N);
            
        end
        
        function u_next = feedforward(obj,trj,x)
            
            dev = x - trj.s;
            h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            ddev = diff(dev')'/h;
            ddev = ddev(1,:);
            dev = dev(1,2:end);                        
    
            % set learning rate
            a_p = 0.5;
            a_d = 0.0;
            u_next = obj.u_last - a_p * dev + a_d * ddev;
            
        end
        
    end
    
end