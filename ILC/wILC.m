% Motor Primitive based Iterative Learning Control 
%
% updates the weights of the DMP
% where Fw is the lifted weight-to-output matrix
% and L is the learning matrix of the ILC update i.e.
% wnext = wlast - L*error;
%
%

classdef wILC < ILC
    
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
        
        % control inputs are not updated for wILC
        u_last
        % Lifted state matrix F (weight-state) and G (state-output)
        Fw
        G
        % Lifted penalty matrices Q and R
        Ql
        Rl
    end
    
    methods
        
        function obj = wILC(model,trj)
                        
            obj.episode = 0;
            obj.color = 'm';
            obj.name = 'weight-based ILC';
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
            
            obj.Fw = zeros(N*dim_x, N*dim_u);
            obj.G = zeros(N*dim_y, N*dim_x);
            obj.Ql = zeros(N*dim_y, N*dim_y);
            obj.Rl = zeros(N*dim_u, N*dim_u);
            
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
                        obj.Fw(vec_x,vec_u) = mat; % on diagonals only B(:,m)
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
                        obj.Fw(vec_x,vec_u) = mat; % on diagonals only B(:,m)
                    end
                end
                
            end
            
            obj.Fw = obj.G * obj.Fw;
            
        end
        
        function feedforward(obj,dmp,traj,y)
            
            r = traj.s;
            dev = y - r;            
            force = dmp.regression(dev);
            %force = dmp.LWR(dev);
            w_change = force.w;
    
            % set learning rate
            beta = 0.01;
            
            w_last = dmp.FORCE.w;
            w_next = w_last - beta * w_change;
            force.w = w_next;
            dmp.setForcing(force);
            
        end
        
    end
    
end