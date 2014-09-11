% Class implementing a basic Kalman filter

classdef Filter < handle
    
    properties
          
        % Kalman filter covariance of the state
        P
        % state to be estimated
        x
        % A matrix of linear system
        A
        % B matrix transforming inputs to state
        B
        % Observation matrix
        C
        % D matrix from inputs to output
        D
        % process noise covariance
        Omega
        % observation noise covariance
        M       
    end
    
    methods
        
        % initialize variances
        function obj = Filter(model,trj,mats)
            
            dimx = model.dim_x;
            N = trj.N - 1;
            obj.P = eye(dimx*N);
            obj.x = zeros(dimx*N,1);
            obj.A = mats.A;
            obj.B = mats.B;
            obj.C = mats.C;
            obj.D = mats.D;
            
            % initialize covariances
            obj.Omega = mats.Omega;
            obj.M = mats.M;
            
        end
        
        % predict state and covariance
        function predict(obj,u)
            
            obj.x = obj.A * obj.x + obj.B * u;
            obj.P = obj.A * obj.P * obj.A + obj.Omega;
        end
        
        % update Kalman filter variance
        function update(obj,dev,u)
            
            % Innovation (output/residual) covariance
            Theta = obj.C * obj.P * obj.C' + obj.M;
            % Optimal Kalman gain
            K = obj.P * obj.C' * inv(Theta); %#ok
            % update state/variance
            obj.P = obj.P - K * obj.C * obj.P;
            obj.x = obj.x + K * (dev(:) - obj.C*obj.x - obj.D * u);
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end