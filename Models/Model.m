% Model superclass acts as an interface
% for holding specific submodel classes.

classdef (Abstract) Model < handle
    
    % Field necessary for all control algorithms
    properties (Abstract)        
        % parameters structure
        PAR
        % constraints structure
        CON
        % fields necessary for simulation and plotting, noise etc.
        SIM
        % cost structure, not necessarily quadratic
        COST
        % observation matrix
        C
    end
    
    % methods to be implemented
    methods (Abstract)
        
        % generate inputs with Trajectory class as output
        Traj = generateInputs(obj,t,ref)
        % set a nominal model
        [xdot,A,B] = nominal(obj,t,x,u)
        % add a disturbance model to the nominal model
        xdot = actual(obj,t,x,u)
        
    end
    
    % methods that can be implemented here in abstract class
    methods (Access = public)
        
        %% one step simulation along the trajectory
        function next = step(obj,t,prev,u,fun)
            
            h = obj.SIM.h;
            dim = obj.SIM.dimx/2;
            
            if obj.SIM.discrete
                next = fun(t,prev,u);            
            elseif strcmpi(obj.SIM.int,'Euler')
                delta = h * fun(t,prev,u);
                next = prev + delta;          
            elseif strcmpi(obj.SIM.int,'Symplectic Euler')
                delta = h * fun(t,prev,u);
                % compute the velocities
                next_vel = prev(dim+1:end) + delta(dim+1:end);
                % compute positions with those velocities
                next_pos = prev(1:dim) + h * next_vel;
                next = [next_pos; next_vel];
            elseif strcmpi(obj.SIM.int,'RK4')
                % get trajectory of states
                % using classical Runge-Kutta method (RK4)  
                k1 = h * fun(t,prev,u);
                x_k1 = prev + k1/2;
                k2 = h * fun(t,x_k1,u);
                x_k2 = prev + k2/2;
                k3 = h * fun(t,x_k2,u);
                x_k3 = prev + k3;
                k4 = h * fun(t,x_k3,u);
                next = prev + (k1 + 2*k2 + 2*k3 + k4)/6;
            else
                error('Method not implemented. Quitting...');
            end

        end
        
        %% Prediction functions
        
        % predict one step using nominal model
        function x_pre = predict(obj,t,x,u)         
            x_pre = step(obj,t,x,u,@obj.nominal);
        end
        
        % predict next states using function
        % useful to calculate nominal model prediction error
        function x_pre = predict_full(obj,t,x,us)
            N = length(t)-1;
            x_pre = zeros(size(x,1),N+1);
            x_pre(:,1) = x(:,1);
            for i = 1:N
                x_pre(:,i+1) = step(obj,t(i),x(:,i),us(:,i),...
                                    @obj.nominal);
            end            
        end
        
        %% Rollout methods here
        
        % useful to propagate feedback law
        % Feedback has to be in state-feedback form!
        function [y,u] = observeWithFeedback(obj,traj,x0)
            fun = @(t,x,u) obj.actual(t,x,u);
            t = traj.t;
            K = traj.K;
            uff = traj.unom;
            N = length(t)-1;
            x = zeros(length(x0),N+1);
            u = zeros(size(K,1),N);
            x(:,1) = x0;
            y(:,1) = obj.C * x(:,1);
            for i = 1:N
                u(:,i) = K(:,:,i)*x(:,i) + uff(:,i);
                x(:,i+1) = step(obj,t(i),x(:,i),u(:,i),fun);
                y(:,i+1) = obj.C * x(:,i+1);
                % no constraint checking
            end
        end
        
        % useful to propagate tILC
        % LQR-calculated K has to be in error-feedback form!
        % DMP can be supplied as an additional argument
        function [y,u] = observeWithFeedbackErrorForm(obj,traj,x0,varargin)
            fun = @(t,x,u) obj.actual(t,x,u);
            hTraj = traj.t(2) - traj.t(1);
            traj = traj.upsample(hTraj/obj.SIM.h);
            t = traj.t;                        
            N = length(t)-1;
            K = traj.K;
            uff = traj.unom;
            
            % process and measurement noise same for now
            sigma = obj.SIM.eps ^(1/2);

            % in case dmp is supplied
            if nargin == 4
                dmps = varargin{:};
                for j = 1:length(dmps)
                    [~,Q] = dmps(j).evolve(N+1);
                    q(j,:) = Q(1,:);
                end
                Cbar = obj.C;
                sbar = Cbar'*((Cbar*Cbar')\q);
            elseif size(traj.s,1) == length(x0)
                sbar = traj.s;
            else
                sbar = traj.projectBack(obj.C);
            end

            x = zeros(length(x0),N+1);
            u = zeros(size(K,1),N);
            x(:,1) = x0;
            y(:,1) = obj.C * x(:,1) + obj.C * sigma * randn(length(x0),1);
            for i = 1:N
                x(:,i) = x(:,i) + obj.C * sigma * randn(length(x0),1);
                u(:,i) = K(:,:,i)*(x(:,i)-sbar(:,i)) + uff(:,i);
                x(:,i+1) = step(obj,t(i),x(:,i),u(:,i),fun);
                y(:,i+1) = obj.C * x(:,i+1);
                % no constraint checking
            end
        end
        
        % useful to propagate one full iteration of
        % ILC input sequence
        function x_next = evolve(obj,t,x0,us)
            %fun = @(t,x,u) obj.actual(t,x,u);
            x_next = simulate(obj,t,x0,us,@obj.actual);
        end
        
        % simulates whole trajectory
        function x = simulate(obj,t,x0,us,fun)
            N = length(t)-1;
            h = t(2)-t(1);
            x = zeros(length(x0),N+1);
            x(:,1) = x0;
            for i = 1:N
                x(:,i+1) = step(obj,t(i),x(:,i),us(:,i),fun);
                % no constraint checking
            end
        end
        
        % add observation noise to state evolution
        % lifted vector representation
        function y = observe(obj,t,x0,us)
            
            N = length(t);
            dimx = obj.SIM.dimx;
            dimy = size(obj.C,1);
            eps = obj.SIM.eps; % covariance of process x measurement noise
            M = eps * eye(N*dimx);
            Cs = cell(1,N);
            [Cs{:}] = deal(obj.C);
            C = blkdiag(Cs{:});
            
            xact = evolve(obj,t,x0,us);
            % vectorize x_iter into N*dim dimensions
            x_vec = xact(:);
            % add observation noise with covariance M
            if eps ~= 0
                x_vec = x_vec + chol(M)*randn(length(x_vec),1);
            end
            % measure only some of the states
            y_vec = C*x_vec;
            % arrange back to normal form
            y = reshape(y_vec,dimy,N);
            
        end
        
        %% Creates a dmp trajectory
        % numbf: number of basis functions to use
        function [dmp,s] = dmpTrajectory(obj,t,numbf,goal,yin,ref)
                                    
            h = obj.SIM.h;
            N = length(t);
            
            % robot state space is of the form [q1,...,qn,qd1,...,qdn];
            if isa(obj,'Robot')
                dim = size(ref,1)/2;
                s = zeros(2*dim,N);
            else
            % for linear models we generate dmps for all the states
                dim = size(ref,1);
                s = zeros(dim,N);
            end
            
            % make a DMP that smoothens x_des
            pat = 'd';
            % scaling
            tau = 1;
            % ensure that phase decays
            exponent = 1;
            ax = exponent/(tau*t(end));            
            can = CAN(h,ax,tau,numbf,t(end),pat);
            
            % create different DMPs, one for each dimension
            alpha = 25;
            beta = 25/4;

            for i = 1:dim
                y0 = [yin(i,:)';0];
                % create the dmp trajectory
                dmp(i) = DDMP(can,alpha,beta,goal(i),y0);
                % learn the weights using regression
                dmp(i).updateWeights(ref(i,:));
                % evolve the DMP
                [x,si] = dmp(i).evolve(N);
                if isa(obj,'Robot')
                    s(i,:) = si(1,:);
                    s(i+dim,:) = si(2,:);
                else
                    s(i,:) = si(1,:);
                end
                
                %figure;
                %plot(t,ref(i,:),'-',t,s(i,:),'-.',t,x);
                %legend('reference trajectory','state y','phase');
                %title('DMP trajectory');
            end         

        end
        
        %% linearizes the nominal dynamics around the trajectory
        function [A,B] = linearize(obj,trj)
            
            N = trj.N - 1; 
            t = trj.t;
            h = trj.t(2) - trj.t(1);
            s = trj.projectBack(obj.C); % get full state space repr.
            
            % if learning takes place in cartesian space then run inverse
            % kinematics
            if isa(obj,'Robot')
                if ~obj.flag_jspace
                    % in this case run inverse kinematics
                    s = obj.invKinematics(s);
                    % add joint velocities
                    sd = diff(s')'/ h; sd(:,end+1) = sd(:,end);
                    s = [s;sd];
                end
            end
            
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            %h = obj.SIM.h;
            unom = trj.unom;
            A = zeros(dimx,dimx,N);
            B = zeros(dimx,dimu,N);
            for i = 1:N
                [~,A(:,:,i), B(:,:,i)] = obj.nominal(t(i),s(:,i),...
                                                 unom(:,i),true);
                % get discrete approximation from jacobian
                % crude approximation
                %A(:,:,i) = eye(dimx,dimx) + obj.h * A(:,:,i);
                %B(:,:,i) = obj.h * B(:,:,i);
                % exact matrix calculation 
                Mat = [A(:,:,i), B(:,:,i); zeros(dimu, dimx + dimu)];
                MD = expm(h * Mat);
                A(:,:,i) = MD(1:dimx,1:dimx);
                B(:,:,i) = MD(1:dimx,dimx+1:end);
            end
        end
        
        %% Plotting methods
        
        % plot the control inputs
        function plot_inputs(obj,trj)
            
            if ~isempty(trj.PERF)
                u = trj.PERF(end).u;
                t = trj.t(1:end-1);
            else
                u = trj.unom;
                t = trj.t(1:end-1);
            end
            
            num_inp = size(u,1);
            figure;
            for i = 1:num_inp
                subplot(num_inp,1,i);
                plot(t,u(i,:));
                title(strcat(num2str(i),'. control input'));
                if ~isa(obj,'BarrettWAM')
                    xlabel('Time (s)');
                    ylabel(strcat('Control input u',num2str(i)));
                end
            end
            
        end
        
        % plot the desired system states and outputs
        function plot_outputs(obj,trj)
 
            t = trj.t;
            y_des = trj.s;
            if ~isempty(trj.PERF)
                y = trj.PERF(end).y;                
                name = trj.PERF(end).name;
            else
                error('Performance not added!');
            end
            num_out = size(y,1);
            figure;
            
            if isa(obj,'BarrettWAM') % too many joints       
                for i = 1:num_out/2
                    subplot(num_out/2,2,2*i-1);
                    plot(t,y(2*i-1,:),'.-',t,y_des(2*i-1,:),'-');
                    subplot(num_out/2,2,2*i);
                    plot(t,y(2*i,:),'.-',t,y_des(2*i,:),'-');
                    %legend(name,'Reference');
                    %title(strcat(num2str(i),'. state'));
                    %xlabel('Time (s)');
                    %ylabel(strcat('State x',num2str(i)));
                end
            else % other models            
                for i = 1:num_out
                    subplot(num_out,1,i);
                    plot(t,y(i,:),'.-',t,y_des(i,:),'-');
                    legend(name,'Reference');
                    title(strcat(num2str(i),'. state'));
                    xlabel('Time (s)');
                    ylabel(strcat('State x',num2str(i)));
                end
            
            end
            
        end
        
    end
end