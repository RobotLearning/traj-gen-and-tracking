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
        
        % inputs are the weights for wILC
        inp_last
        % Lifted state matrix F (weight-state) and G (state-output)
        Fw
        G
        % Lifted penalty matrix Q
        Ql
        % Weight penalty matrix R
        Rw
    end
    
    methods
        
        function obj = wILC(model,dmp,trj)
                        
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
            
            % number of basis functions
            m = length(dmp.FORCE.c);
            
            obj.Fw = zeros(4*N*dim_x, m);
            obj.G = zeros(N*dim_y, N*dim_x);
            obj.Ql = zeros(N*dim_y, N*dim_y);
            obj.Rw = zeros(m, m);
            
            % fill the F matrix
            %obj.lift(model,dmp,trj);
            
        end
        
        % get the lifted vector representation 
        % around the trajectory
        function lift(obj,model,dmp,trj)
            
            N = trj.N - 1;
            m = length(dmp.FORCE.c);
            dim_x = model.SIM.dimx;
            dim_u = model.SIM.dimu;
            C = [eye(dim_x), -eye(dim_x), zeros(dim_x,2*dim_x)];
            
            % deal C matrix to G
            obj.G = cell(1,N);
            obj.Ql = cell(1,N);
            obj.Rw = cell(1,m);
            % TODO: this could also be time varying!
            [obj.G{:}] = deal(model.C);
            [obj.Ql{:}] = deal(model.COST.Q);
            [obj.Rw{:}] = deal(model.COST.R);
            obj.G = blkdiag(obj.G{:});
            obj.Ql = blkdiag(obj.Ql{:});
            obj.Rw = blkdiag(obj.Rw{:});
            
            % transition matrices
            obj.Apsi = zeros(4*N*dim_x,4*N*dim_x);
            obj.Bpsi = zeros(4*N*dim_x,m);
            
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
            
            % appears because of static update of the weights
            mat = repmat(eye(m),N,1);
            obj.Fw = obj.G * obj.Fw * mat;
            
        end
        
        function feedforward(obj,dmp,traj,y)
            
            r = traj.s;
            dev = y - r;
            %dev = zeros(1,length(dev));
            force = dmp.updateWeights(dev);
            w_change = force.w;
    
            % set learning rate
            beta = 0.01;
            
            w_last = dmp.FORCE.w;
            w_next = w_last + beta * w_change;
            force.w = w_next;
            dmp.setForcing(force);
            
        end
        
    end
    
end