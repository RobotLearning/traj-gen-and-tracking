% Abstract Iterative Learning Control (ILC) superclass 
% which all ILC implementations must inherit from
% for holding performances, plotting, etc...

classdef (Abstract) ILC < handle
    
    properties (Abstract)
        
        % number of total episodes so far
        episode
        % color of particular controller
        color
        % name of the particular controller
        name
        % RMS error incurred along the trajectory
        error
        % final cost
        finalCost
        % downsampling to speed things up
        downsample        
        % ILC's Last input sequence
        inp_last
        
    end
    
    methods (Abstract)
        
        % apply feedforward control (sequence of control signals)
        feedforward(i,trj,model,x)
        
    end
    
    methods
        
        % get sse costs for the controller
        function record(obj,inp_applied,err,finCost)            
            obj.episode = obj.episode + 1;
            obj.inp_last = inp_applied;
            obj.error(obj.episode) = err;          
            obj.finalCost(obj.episode) = finCost;
        end
        
        % downsampling function to speed things up
        %function downsample(obj,inp_applied)
        %    %TODO
        %end
        
    end
    
    
    
end