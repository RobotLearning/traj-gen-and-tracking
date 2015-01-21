% Model based Iterative Learning Control 
%
% where F is the lifted matrix of the plant dynamics
% and L is the learning matrix of the ILC update i.e.
% unext = ulast - L*error;
%
% Ideally we need L = inv(F), a model-based ILC update rule

classdef mILC < ILC
    
    % fields common to all ILCs (defined in abstract ILC class)
    properties
         
        % number of total episodes so far
        episode
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        error
        % flag for learning with feedback? [then model changes]
        learn_fb
        
        % ILC's Last input sequence
        u_last
        % Lifted state matrix F (input-state) and G (state-output)
        F
        G
        % Lifted penalty matrices Q and R
        Ql
        Rl
    end
    
    methods
        
        % varargin is for learning with feedback, better way?
        function obj = mILC(model,trj,varargin)
                        
            obj.episode = 0;
            obj.color = 'm';
            obj.name = 'Model-based ILC';
            obj.error = 0;
            
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            dim_y = model.SIM.dimy;
            
            N = trj.N - 1;
            if ~isempty(trj.unom)
                obj.u_last = trj.unom(:,1:N);
            else
                warning('Using last performance results');
                obj.u_last = trj.PERF(end).u;
            end
            
            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = zeros(N*dim_y, N*dim_x);
            obj.Ql = zeros(N*dim_y, N*dim_y);
            obj.Rl = zeros(N*dim_u, N*dim_u);
            
            obj.learn_fb = false;
            if nargin == 3
                obj.learn_fb = true;
            end
            
            % fill the F matrix
            obj.lift(model,trj);
            
        end
        
        % get the lifted vector representation 
        % around the trajectory
        function lift(obj,model,trj)
            
            N = trj.N - 1;
            
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            
            % deal C matrix to G
            obj.G = cell(1,N);
            obj.Ql = cell(1,N);
            obj.Rl = cell(1,N);
            % TODO: this could also be time varying!
            [obj.G{:}] = deal(model.C);
            [obj.Ql{:}] = deal(model.COST.Q);
            [obj.Rl{:}] = deal(model.COST.R);
            obj.G = blkdiag(obj.G{:});
            obj.Ql = blkdiag(obj.Ql{:});
            obj.Rl = blkdiag(obj.Rl{:});
            
            if isa(model,'Linear')
            
                Ad = model.Ad;
                Bd = model.Bd;
                if ~obj.learn_fb
                    % construct lifted domain matrix F
                    % TODO: this can be computed much more efficiently
                    % for linear systems
                    for i = 1:N
                        for j = 1:i
                            vec_x = (i-1)*dim_x + 1:i*dim_x;
                            vec_u = (j-1)*dim_u + 1:j*dim_u;
                            mat = Bd;
                            for k = j+1:i
                                mat = Ad * mat;
                            end
                            obj.F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                        end
                    end
                else % learn with feedback
                    K = trj.K;
                    % construct lifted domain matrix F
                    % TODO: this can be computed much more efficiently
                    % for linear systems
                    for i = 1:N
                        for j = 1:i
                            vec_x = (i-1)*dim_x + 1:i*dim_x;
                            vec_u = (j-1)*dim_u + 1:j*dim_u;
                            mat = Bd;
                            for k = j+1:i
                                mat = (Ad + Bd * K(:,:,k)) * mat;
                            end
                            obj.F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                        end
                    end
                end
            else
                % get linear time variant matrices around trajectory
                [Ad,Bd] = model.linearize(trj);
                
                if ~obj.learn_fb
                    % construct lifted domain matrix F
                    for i = 1:N
                        for j = 1:i
                            vec_x = (i-1)*dim_x + 1:i*dim_x;
                            vec_u = (j-1)*dim_u + 1:j*dim_u;
                            mat = Bd(:,:,j);
                            for k = j+1:i
                                mat = Ad(:,:,k) * mat;
                            end
                            obj.F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                        end
                    end
                else % learn with feedback
                    K = trj.K;
                    % construct lifted domain matrix F
                    for i = 1:N
                        for j = 1:i
                            vec_x = (i-1)*dim_x + 1:i*dim_x;
                            vec_u = (j-1)*dim_u + 1:j*dim_u;
                            mat = Bd(:,:,j);
                            for k = j+1:i
                                mat = (Ad(:,:,k) + Bd(:,:,k)*K(:,:,k)) * mat;
                            end
                            obj.F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                        end
                    end
                end
            end
            
            obj.F = obj.G * obj.F;
            
        end
        
        function u = feedforward(obj,trj,y)
            
            dimu = size(obj.u_last,1);
            N = size(obj.u_last,2);
            
            dev = y - trj.s;
            h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            % ddev = diff(dev')'/h;
            % ddev = ddev(1,:);
            dev = dev(:,2:end);                        
    
            % set learning rate
            beta = 10;
            
            % gradient descent
            %u = obj.u_last(:) - beta * obj.F' * obj.Ql * dev(:);
            % model inversion based Newton-Raphson update
            %u = obj.u_last(:) - obj.F \ dev(:);
            % more stable inverse based Newton-Raphson update
            % computes very high inverses though
            u = obj.u_last(:) - pinv(obj.F) * dev(:);
            % LM-type update
            %Mat = (obj.F' * obj.Ql * obj.F + obj.Rl) \ (obj.F' * obj.Ql);
            %u = obj.u_last(:) - Mat * dev(:);
            
            % revert from lifted vector from back to normal form
            u = reshape(u,dimu,N);
            
        end
        
    end
    
end