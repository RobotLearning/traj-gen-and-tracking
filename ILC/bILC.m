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
        inp_last
    end
    
    methods
        
        %% Constructor for basic Arimoto-type ILC modifying u_ff
        function obj = bILC(traj)
                        
            obj.episode = 0;
            obj.color = 'b';
            obj.name = 'Arimoto type ILC';
            obj.error = 0;
            
            N = traj.N - 1;
            if ~isempty(traj.unom)
                obj.inp_last = traj.unom(:,1:N);
            else
                warning('Using last performance results');
                obj.inp_last = traj.PERF(end).u;
            end
            
        end
        
        %% Main ILC function
        function u_next = feedforward(obj,traj,y)
            
            dev = y - traj.s;
            %h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            ddev = diff(dev')';
            dev = dev(:,2:end);                        
    
            % set learning rate
            a_p = 0.1;
            a_d = 0.2;
            u_next = obj.inp_last - a_p * dev - a_d * ddev;
            
        end
        
        % version for robot
        function u_next = feedforwardRobot(obj,traj,y)
            
            dev = y - traj.s;
            %h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            dimu = size(dev,1)/2;
            e = dev(1:dimu,2:end);                        
            edot = dev(dimu+1:end,2:end);
    
            % set learning rate
            a_p = 0.05;
            a_d = 0.05;
            u_next = obj.inp_last - a_p * e - a_d * edot;
        end
        
    end
    
end