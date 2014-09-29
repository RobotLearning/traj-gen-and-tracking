% Discrete Dynamic Movement Primitive

classdef discreteDMP < DMP
    
    properties
        
    end
    
    methods
        
        function obj = discreteDMP()
            
        end
        
        % basis functions are unscaled gaussians
        function out = basis(x,h,c)
        out = exp(-h * (x - c)^2);
        end
        
        % tests the dmp
        function out = test()
            
        end
        
    end
    
end