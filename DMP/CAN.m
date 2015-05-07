% Canonical system which is being shared by multiple DMPs

classdef CAN < handle
    
    properties
        
        % phase variable which synchronizes DMPs
        x
        % time step constant
        a_x
        % temporal scaling factor
        tau
        % boolean pattern variable - discrete or rhythmic
        pattern
        % time step
        dt
        % number of basis functions
        nbf
        % centers of the forcing function rbfs
        c
        % heights of the forcing function rbfs
        h
    end
    
    methods (Access = public)
        
        %% Constructor for the common phase of DMPs
        function obj = CAN(h, ax, tau, bfs, ttotal, pat)
        
            assert(strcmp(pat,'d') || strcmp(pat, 'r'), ...
                   'Please input d for discrete and r for rhythmic');
            obj.pattern = pat;          
            obj.dt = h;
            obj.a_x = ax;
            obj.tau = tau;
            
            % initialize the forcing function parameters
            obj.nbf = bfs;
            obj.c = zeros(1,bfs);
            obj.h = zeros(1,bfs);
            
            t_total_guess = ttotal; % guess running time
            t = (t_total_guess/(bfs-1)) * ((1:bfs)-1);
            for i = 1:bfs
                obj.c(i) = exp(-ax * t(i));
            end
            
            for i = 1:bfs-1
                obj.h(i) = 0.5 / ((0.65 * (obj.c(i+1)-obj.c(i))) ^ 2);
            end
            obj.h(end) = obj.h(end-1);
            
            % call reset to make sure we start from zero
            obj.reset();
        end
        
        %% Methods modifying the phase
        
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
        function x_tr = evolve(obj,N)
            
            % make sure phase is reset
            obj.reset();
            x_tr = zeros(1,N);
            for i = 1:N
                x_tr(1,i) = obj.x;
                obj.step(1);
            end  
        end
        
    end

end