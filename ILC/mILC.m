% Monotonic Iterative Learning Control 
% that tries to satisfy the monotonic update rule:
%
% max. singular value s(I - LF) < 1
%
% where F is the lifted matrix of the plant dynamics
% and L is the learning matrix of the ILC update i.e.
% unext = ulast + L*error;
%
% Ideally we would need L = inv(F), a model-based ILC update rule

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
        
        % ILC's Last input sequence
        u_last
        % Lifted state matrix F (input-state) and G (state-output)
        F
        G
    end
    
    methods
        
        function obj = mILC(model,trj)
                        
            obj.episode = 0;
            obj.color = 'm';
            obj.name = 'Monotonic ILC';
            obj.error = 0;
            
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            dim_y = model.SIM.dimy;
            
            N = trj.N - 1;
            obj.u_last = trj.unom(:,1:N);
            
            obj.F = zeros(N*dim_x, N*dim_u);
            obj.G = zeros(N*dim_y, N*dim_x);
            
            % fill the F matrix
            obj = obj.lift(model,trj);
            
        end
        
        % get the lifted vector representation 
        % around the trajectory
        function obj = lift(obj,model,trj)
            
            N = trj.N - 1;
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            
            % deal C matrix to G
            obj.G = cell(1,N);
            % TODO: this could also be time varying!
            [obj.G{:}] = deal(model.C);
            obj.G = blkdiag(obj.G{:});
            
            if isa(model,'Linear')
                
                Ad = model.Ad;
                Bd = model.Bd;
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

                
            else
                % get linear time variant matrices around trajectory
                [Ad,Bd] = model.linearize(trj);
                
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
                
            end
            
            obj.F = obj.G * obj.F;
            
        end
        
        function u_next = feedforward(obj,trj,y)
            
            dev = y - trj.s;
            h = trj.t(2) - trj.t(1);
            % get rid of x0 in dev
            % ddev = diff(dev')'/h;
            % ddev = ddev(1,:);
            dev = dev(1,2:end);                        
    
            % set learning rate
            % a_p = 0.5;
            % a_d = 0.2;
            % u_next = obj.u_last - a_p * dev - a_d * ddev;
            
            % model based update
            u_next = obj.u_last(:) - obj.F \ dev(:);
            u_next = u_next';
            
        end
        
    end
    
end