% Rhythmic Dynamic Movement Primitive

classdef rhythmicDMP < DMP
    
    properties
        
    end
    
    methods
        
        % basis functions are unscaled von-mises functions 
        % i.e. periodic gaussian
        function out = vonMises(phi,h,c)
        out = exp(h * (cos(phi - c) - 1));
        end
    end
    
end