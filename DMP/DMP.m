% Dynamic motor primitive superclass
% rhythmic and discrete DMPs inherit from this class

classdef (Abstract) DMP < handle
    
    properties (Abstract)
        
        % time constants
        alpha_g, beta_g
        % goal state
        g
        % temporal scaling factor
        tau
        % time evolution as a [t0,tf] vector
        tspan
        % initial y, yd (or z) values
        y0, yd0
    end
    
    methods (Abstract)
        
        % canonical system - differs for rhytmic and discrete DMP
        canonical()        
        
        % discrete and rhythmic subclasses define their own basis fnc
        basis(phi,h,c)
        
        % make a step - both for feedback/feedforward simulation
        step()
        
        % learn weights using this function
        % basis functions are fixed
        learnWeightsFixed()
        
        % phase of the canonical system is set to 0
        resetPhase()
        
        % subclasses must implement a -simple- test scenario
        test()
        
    end    
    
end