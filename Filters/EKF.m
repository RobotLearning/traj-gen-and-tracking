% Class implementing an Extended Kalman filter
% Assuming noise covariances are time invariant

classdef EKF < handle
    
    properties
          
        % Kalman filter covariance of the state
        P
        % state to be estimated
        x
        % function handle for state evolution model
        f
        % function handle for observation model
        g
        % linearization of the function f around estimated state
        A
        % linearization of the function g around estimated state
        C
        % process noise covariance
        Q
        % observation noise covariance
        R
    end
    
    methods
        
        %% Initialize variances
        function obj = EKF(dim,funState,funObserve,mats)
            
            dimx = dim;
            obj.P = eye(dimx);
            obj.x = zeros(dimx,1);
            obj.f = funState;
            obj.g = funObserve;
            
            % initialize covariances
            obj.Q = mats.O;
            obj.R = mats.M;
            
            % initialize linearized matrices
            obj.A = [];
            obj.C = [];
            
        end
        
        %% linearize the dynamics f and g with 'TWO-SIDED-SECANT'
        % make sure you run this before predict and update
        function linearize(obj)
            
            % TWO-SIDED-SECANT
            % Take 2nd differences with h small
            h = 1e-3;
            n = length(obj.x); 
            xhPlus = repmat(obj.x,1,2*n) + h * eye(2*n);
            xhMinus = repmat(obj.x,1,2*n) - h * eye(2*n);
            fPlus = zeros(2*n); 
            gPlus = zeros(2*n);
            fMinus = zeros(2*n);
            gMinus = zeros(2*n);
            for i = 1:2*n
                fPlus(:,i) = obj.f(xhPlus(:,i),u);
                fMinus(:,i) = obj.f(xhMinus(:,i),u);
                gPlus(:,i) = obj.g(xhPlus(:,i),u);
                gMinus(:,i) = obj.g(xhMinus(:,i),u);
            end
            fder = (fPlus - fMinus) / (2*h);
            dgdx = (gPlus - gMinus) / (2*h);
            dfdx = [zeros(n), eye(n); fder(n+1:2*n,:)];
            
            obj.A = dfdx;
            obj.C = dgdx;
        end
        
        %% Initialize state
        function initState(obj,x0,P0)
            obj.x = x0;
            obj.P = P0;
        end
        
        %% predict state and covariance
        function predict(obj,u)
            
            obj.x = obj.f(obj.x,u);
            obj.P = obj.A * obj.P * obj.A + obj.Q;
        end
        
        %% update Kalman filter variance
        function update(obj,y,u)
            
            % Innovation sequence
            Inno = y(:) - obj.g(obj.x,u);
            % Innovation (output/residual) covariance
            Theta = obj.C * obj.P * obj.C' + obj.R;
            % Optimal Kalman gain
            K = obj.P * obj.C' * inv(Theta); %#ok
            % update state/variance
            obj.P = obj.P - K * obj.C * obj.P;
            obj.x = obj.x + K * Inno;
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end