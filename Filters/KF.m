% Class implementing a basic Kalman filter
% assuming an LTI system

classdef KF < handle
    
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
        Q
        % observation noise covariance
        R
    end
    
    methods
        
        %% Initialize variances
        function obj = KF(dim,mats)
            
            dimx = dim;
            obj.P = eye(dimx);
            obj.x = zeros(dimx,1);
            obj.A = mats.A;
            obj.B = mats.B;
            obj.C = mats.C;
            obj.D = mats.D;
            
            % initialize covariances
            obj.Q = mats.O;
            obj.R = mats.M;
            
        end
        
        %% Initialize state
        function initState(obj,x0,P0)
            obj.x = x0;
            obj.P = P0;
        end
        
        %% predict state and covariance
        function predict(obj,u)
            
            obj.x = obj.A * obj.x + obj.B * u;
            obj.P = obj.A * obj.P * obj.A' + obj.Q;
        end
        
        %% update Kalman filter variance
        function update(obj,y,u)
            
            % Innovation sequence
            Inno = y(:) - obj.C*obj.x - obj.D*u;
            % Innovation (output/residual) covariance
            Theta = obj.C * obj.P * obj.C' + obj.R;
            % Optimal Kalman gain
            K = obj.P * obj.C' * inv(Theta); %#ok
            % update state/variance
            obj.P = obj.P - K * obj.C * obj.P;
            obj.x = obj.x + K * Inno;
        end
        
        %% Kalman smoother (batch mode)
        % class must be initialized and initState must be run before
        function [X,V] = smooth(obj,Y,U)
            
            N = size(Y,2);
            X = zeros(length(obj.x),N);
            X_pred = zeros(length(obj.x),N-1);
            V = zeros(length(obj.x),length(obj.x),N);
            V_pred = zeros(length(obj.x),length(obj.x),N-1);
            % forward pass
            for i = 1:N-1
                obj.update(Y(:,i),U(:,i));
                X(:,i) = obj.x;
                V(:,:,i) = obj.P;
                obj.predict(U(:,i));
                X_pred(:,i) = obj.x;
                V_pred(:,:,i) = obj.P;
            end
            obj.update(Y(:,N),0);
            X(:,N) = obj.x;
            V(:,:,N) = obj.P;
            % backward pass
            for i = N-1:-1:1
                % Rauch recursion 
                H = V(:,:,i) * obj.A' * inv(V_pred(:,:,i));
                X(:,i) = X(:,i) + H * (X(:,i+1) - X_pred(:,i));
                V(:,:,i) = V(:,:,i) + H * (V(:,:,i+1) - V_pred(:,:,i)) * H';
            end
        end
        
    end
    
end