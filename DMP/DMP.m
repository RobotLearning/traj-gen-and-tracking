% Dynamic motor primitive superclass
% rhythmic and discrete DMPs inherit from this class

classdef (Abstract) DMP < handle
    
    properties (Abstract)
        
        % canonical system
        can
        % time constants
        alpha_g, beta_g
        % goal state
        goal
        % y, yd (or z), ydd values
        Y
        % forcing structure has weights w, widths h, and centers c
        FOR
    end
    
    methods (Abstract)

        % discrete and rhythmic subclasses define their own basis fnc
        basis(phi,h,c)
        
        % forcing function
        forcing(obj)
        
        % make a step - both for feedback/feedforward simulation
        step(obj)
        
        % feedforward simulation of DMPs
        evolve(obj)
        
        % phase of the canonical system is set to 0
        % as well as clearing Y values
        resetStates(obj)
        
        % set the initial state of the DMP
        setInitState(obj,Y0)
        
        % set the goal state from observational data
        setGoal(obj,path)
        
        % set the forcing function after regression for instance
        setForcing(obj,force)
        
    end    
    
end