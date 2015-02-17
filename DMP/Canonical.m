% Canonical system which is being shared by multiple DMPs

classdef Canonical < handle
    
    properties
        
        % phase variable which synchronizes DMPs
        x
        % time step constant
        a_x
        % temporal scaling factor
        tau
        % boolean pattern variable - discrete or rhythmic
        pattern
        % number of time steps
        N
        % time step
        dt
    end
    
    methods (Access = public)
       
        % TODO: provide default values for ax and tau
        function obj = Canonical(h, ax, tau, num_steps, pat)
        
            assert(strcmp(pat,'d') || strcmp(pat, 'r'), ...
                   'Please input d for discrete and r for rhythmic');
            obj.pattern = pat;
            obj.N = num_steps;            
            obj.dt = h;
            obj.a_x = ax;
            obj.tau = tau;
            
            % call reset to make sure we start from zero
            obj.reset();
        end
        
        function reset(obj)
            
            obj.x = 1;
        end
        
        % err is the error coupling term due to feedback error
        % slows down the canonical system 
        % can be used for varying dt as well
        function step(obj,err)
           
            if strcmp(obj.pattern,'d')
                obj.x = obj.x - (obj.tau * obj.a_x * obj.x * err) * obj.dt;
            else
                % assuming rhythmic dmp
                obj.x = obj.x + (err * obj.tau) * obj.dt;
            end
            
        end
        
        % evolve is the feedforward rollout function
        % TODO: apply bsxfun or arrayfun
        function x_tr = evolve(obj)
            
            % make sure phase is reset
            obj.reset();
            x_tr = zeros(1,obj.N);
            for i = 1:obj.N
                x_tr(1,i) = obj.x;
                obj.step(1);
            end  
        end
        
    end

end