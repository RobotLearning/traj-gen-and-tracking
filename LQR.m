% Class for calculating linear feedback law
% TODO: include Kalman filter for LQG
%
% Discrete time finite horizon linear quadratic regulator
% lqr function gives the optimal linear feedback law
%
% Inputs:
% Q,R - Weight matrices
% A,B - Matrices of the differential equation
% Functions - Infinite horizon or finite horizon, as well as time
% varying/invariant
%
% Outputs:
% K - Feedback gain
%
% TODO: extend to include feedforward term
%       Qf can be also input 

classdef LQR 

    properties   
        % State penalty matrix
        Q
        % Input penalty matrix
        R
        % Time horizon
        N
        % discretization to get Ad, Bd
        h
        % A and B matrices
        A, B
        % discrete A and B matrices
        Ad, Bd
        % observation matrix
        C
        % final state penalty matrix
        Qf
    end

    methods
        
        % constructor for convenience
        % TODO: put asserts
        function obj = LQR(Q,R,Qf,A,B,C,N,h)
            
            % construct Q2 - matrix of nonoutput states
            dimy = size(Q,1);
            dimx = size(B,1);
            % all nonoutput states should go to zero
            % TODO: this can be modified
            Q2 = diag([zeros(1,dimy), ones(1,dimx-dimy)]);
            Cbar = eye(dimx) - C'*((C*C')\C);
            obj.Q = Cbar'*Q2*Cbar + C'*Q*C;
            obj.R = R;
            obj.Qf = Cbar'*Q2*Cbar + C'*Qf*C;
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.N = N;
            obj.h = h;
            
            % get Ad, Bd matrices
            obj = obj.discretizeMatrices();
            
        end
        
        function obj = discretizeMatrices(obj)
            
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            % trick to get discrete time versions
            Mat = [obj.A, obj.B; zeros(dimu, dimx + dimu)];
            MD = expm(obj.h * Mat);
            obj.Ad = MD(1:dimx,1:dimx);
            obj.Bd = MD(1:dimx,dimx+1:end);
        end
        
        % TODO: replace with output controllability
        % make sure the system is controllable/reachable
        % otherwise give an error
        function assertControllability(obj)
            
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            % construct controllability Kalman matrix
            K = zeros(dimx,dimx*dimu);
            for i = 0:dimx-1
                K(:,(dimx*i+1):(dimx*i+dimu)) = (obj.Ad^i)*obj.Bd;
            end
            assert(rank(K) == dimx, 'System is not controllable!');
        end
        
        % TODO:
        function assertDetectability(obj)
            
        end
        
        % assumed to be LTI
        function K = computeInfHorizonLTI(obj)
            
            % inf horizon
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            Q = obj.Q;
            R = obj.R;
            A = obj.Ad;
            B = obj.Bd;
            
            Pinf = Q;
            % iterate till convergence
            delta = 1e-3;
            diff = 1;
            i = 0;
            while diff > delta
                Ppre = Pinf;
                Pinf = dynamicRiccati(Pinf,Q,R,A,B);
                matDiff = Ppre - Pinf;
                diff = norm(matDiff,'fro');
                i = i+1;
            end
            fprintf('Iterated dynamic riccati equation %d times\n', i);
            K = -(R + B'*Pinf*B)\(B'*Pinf*A);
        end
        
        % simplest discrete-time case, i.e. no tracking
        function K = computeFinHorizonLTI(obj)
            
            % finite horizon
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            Q = obj.Q;
            R = obj.R;
            A = obj.Ad;
            B = obj.Bd;
            Qf = obj.Qf;
            N = obj.N;

            P = zeros(dimx,dimx,N);
            P(:,:,end) = Qf;

            % calculate optimal K
            K = zeros(dimu,dimx,N);
            K(:,:,end) = -(R + B'*Qf*B)\(B'*Qf*A);
            for i = N-1:-1:1
                P(:,:,i) = dynamicRiccati(P(:,:,i+1),Q,R,A,B);
                K(:,:,i) = -(R + B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
            end
        end
        
        % assumed to be finite horizon (inf. horizon case?)
        function K = computeFinHorizonLTV(obj)
            
            % Arrays A and B hold time-varying matrices
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            Q = obj.Q;
            R = obj.R;
            A = obj.Ad;
            B = obj.Bd;
            Qf = obj.Qf;
            N = obj.N;
            
            P = zeros(dimx,dimx,N);
            P(:,:,end) = Qf;

            % calculate optimal K
            K = zeros(dimu,dimx,obj.N);
            K(:,:,end) = -(R + B(:,:,end)'*Qf*B(:,:,end))\...
                          (B(:,:,end)'*Qf*A(:,:,end));
            for i = N-1:-1:1
                P(:,:,i) = dynamicRiccati(P(:,:,i+1),Q,R,A(:,:,i),B(:,:,i));
                K(:,:,i) = -(R + B(:,:,i)'*P(:,:,i)*B(:,:,i))\...
                            (B(:,:,i)'*P(:,:,i)*A(:,:,i));
            end
        end
        
        % discrete-time finite-horizon trajectory tracking
        % The state-form of LQR tracking is used
        % 
        % Outputs:
        % K and uff calculated using the error form
        function [K,uff] = computeFinHorizonTracking(obj,s)
            
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            Q = obj.Q;
            R = obj.R;
            Qf = obj.Qf;
            A = obj.Ad;
            B = obj.Bd;
            N = obj.N;
            
            P = zeros(dimx,dimx,N);
            P(:,:,end) = Qf;

            % calculate optimal K and uff
            K = zeros(dimu,dimx,N);
            uff = zeros(dimu,N);
            nu = zeros(dimx,N);
            uff(:,end) = 0;
            nu(:,end) = 0;
            K(:,:,end) = -(R + B'*Qf*B)\(B'*Qf*A);
            for i = N-1:-1:1
                P(:,:,i) = dynamicRiccati(P(:,:,i+1),Q,R,A,B);
                K(:,:,i) = -(R + B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
                nu(:,i) = (A + B*K(:,:,i))'*nu(:,i+1) - Q*s(:,i);
                uff(:,i) = -(R + B'*P(:,:,i)*B)\(B'*nu(:,i));
            end
            
        end
        
        % TODO: is this equivalent to error form?
        % discrete-time finite-horizon trajectory tracking
        % The error-form of LQR tracking is used
        % 
        % Outputs:
        % Kbar is K and uff combined in a n+1 x n+1 feedback/feedforward
        % matrix
        function Kbar = computeFinHorizonTrackingStateFeedback(obj,s)
            
            warning('Deprecated function!');
            
            dimu = size(obj.B,2);
            dimx = size(obj.B,1);
            A = obj.Ad;
            B = obj.Bd;
            N = obj.N;
            
            % form the time varying matrices Abar and Bbar
            Abar = zeros(dimx+1,dimx+1,N);
            Bbar = zeros(dimx+1,dimu,N);
            for i = 1:N
                Abar(:,:,i) = [A, A*s(:,i) - s(:,i+1); ...
                               zeros(1,dimx), 1];
                Bbar(:,:,i) = [B; 0];
            end
            
            obj.Ad = Abar;
            obj.Bd = Bbar;
            % modify Q and Qf
            obj.Q = diag([diag(obj.Q);0]);
            obj.Qf = diag([diag(obj.Qf);0]);
            
            Kbar = obj.computeFinHorizonLTV();
            
        end
        
    end
end