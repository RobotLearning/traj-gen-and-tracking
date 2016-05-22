% Class implementing an simple moving average filter
% To estimate position and velocity of a moving target

classdef PolyFilter < handle
    
    properties
          
        % Number of observations to average over
        size
        % Number of observations so far
        numobs
        % Observations to average over
        Y
        % Times of the observations
        T
        % Vandermonde matrix
        V
        % Order of the polynomial
        order
        % Prediction model
        f
        % State (position and velocity)
        x
        % current time
        time
    end
    
    methods
        
        %% Initialize filtersize
        function obj = PolyFilter(dim,order,filter_size,predict_model)
            
            obj.order = order;
            obj.numobs = 0;
            obj.size = filter_size;
            obj.x = zeros(2*dim,1);
            obj.f = predict_model;
            obj.Y = zeros(filter_size,dim);
            obj.V = zeros(filter_size,order+1);
            obj.time = 0.0;
            obj.T = zeros(filter_size,1);
        end
 
        
        %% Initialize state
        function initState(obj,x0,eps)
            assert(length(x0) == 2*size(obj.Y,2),'dim size does not match with initialization!');
            assert(length(x0) == length(obj.x),'x0 size does not match with state!');
            obj.x = x0;
            obj.T = zeros(obj.size,1);
            obj.Y(1,:) = x0(1:size(obj.Y,2));
            obj.numobs = 1;
            obj.time = 0.0;
        end
        
        %% predict state
        % predict dt into the future
        function predict(obj,dt,u)
            
            % integrating the state to get x(t+1)
            obj.x = obj.f(obj.x,dt);
        end
        
        %% update state from observation
        
        % delta_t is the time passed from last observation time
        function update(obj,delta_t,y)
            
            obj.time = obj.time + delta_t;
            if obj.numobs < obj.size
                obj.numobs = obj.numobs + 1;
                obj.Y(obj.numobs,:) = y(:)';
                obj.T(obj.numobs) = obj.time;
                if obj.numobs == obj.size
                    obj.construct_vandermonde();
                    obj.estimate_state();
                end
            elseif obj.numobs == obj.size
                % regress over balls 
                % construct the matrix
                obj.Y = [obj.Y(2:end,:); y(:)'];
                obj.update_vandermonde();
                obj.estimate_state();
            else
                error('this shouldnt happen!');
            end
        end
        
        % construct v-matrix only once
        function construct_vandermonde(obj)
            
            assert(obj.numobs == obj.size, 'numobs should be equal to maxsize');
            for i = 0:1:obj.order
                obj.V(:,i+1) = (obj.T).^i;
            end
        end
        
        % update the vandermonde matrix
        function update_vandermonde(obj)
            
            assert(obj.numobs == obj.size, 'numobs should be equal to maxsize');
            vec = zeros(1,obj.order+1);
            for i = 0:1:obj.order
                vec(i+1) = obj.time^i;
            end
            obj.V = [obj.V(2:end,:);vec];
        end
        
        % Estimate the state from the last 12 observations 
        % after having constructed the vandermonde matrix
        function estimate_state(obj)
             a = obj.V \ obj.Y;
             pos = obj.V(end,:) * a;
             
             vec = zeros(1,obj.order+1);
             for i = 1:obj.order
                 vec(i+1) = factorial(i) * obj.time^(i-1);
             end
             vel = vec * a;
             
             obj.x = [pos(:);vel(:)];
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end