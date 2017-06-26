% Class implementing a Gaussian Process regression model

classdef GP < handle
    
    properties          

        % Hyperparameters
        hp
        % Training data
        x,y
        % Covariance matrix between (training) data points
        cov
        % Kernel function (handle)
        kernel
        
    end
    
    methods (Access = public)
        
        %% Initialize data points, kernel structure and covariance matrix
        function obj = GP(hp,xs,ys)
            
            obj.x = xs;
            obj.y = ys;
            obj.hp = hp;
            obj.setKernel(hp);
            obj.buildCov(xs);
            
        end
        
        % Estimate hyperparameters using maximum likelihood
        function fitHP(obj,hp0)
             
            try
                addpath(genpath('../gpml-matlab-v3.1-2010-09-27'));
            catch
                error('GPML toolbox not found!');
            end
            
            covar = @covSEard;
            lik = @likGauss;
            inf = @infExact;
            mean = [];
            cov_hp_len = size(obj.x,1) + 1; % scale hp added
            hyp.mean = [];
            hyp.cov = [log(hp0.l);log(sqrt(hp0.scale))];
            hyp.lik = log(sqrt(hp0.noise.var));            
            Ncg = 1000; % 100 line search steps
            hyp = minimize(hyp, @gp, Ncg, ...
                              inf, mean, covar, lik, obj.x', obj.y);
                          
            % update hp field after hp estimation
            obj.hp.l = exp(hyp.cov(1:end-1));
            obj.hp.scale = exp(hyp.cov(end))^2;
            obj.hp.noise.var = exp(hyp.lik)^2;
            obj.setKernel(obj.hp);
        end
        
        
        %% predict mean and variance at a test point
        function [mu,s2] = predict(obj,xstar)
            
            % calculate kernel vector k and store
            vec = obj.covTestAndData(xstar);
            % create mu function 
            mu = vec' * (obj.cov \ obj.y);
            % initial covar is 1 + noise
            kxx = obj.kernel(xstar,xstar) + obj.hp.noise.var;
            % subtract information gain
            s2 = kxx - vec' * (obj.cov \ vec);
        end
        
        %% update mean and variance
        function update(obj,xnew,ynew)
            
            obj.x = [obj.x, xnew];
            obj.y = [obj.y; ynew];
            obj.updateCov(xnew);
        end        
    end
    
    methods (Access = private)
        
        %% Kernel related methods
        
        % Kernel called by the ker_matrix, ker_matrix_iter and ker_vector
        % functions.
        % x1 - Vector 1
        % x2 - Vector 2 
        %
        % Sets up a Gaussian or linear kernel function depending on 
        % the type field and hyperparameters
        %
        % out - covariance between x1 and x2
        %
        % IMPORTANT: does not consider noise as part of kernel
        %
        function setKernel(obj,hp)

            L = hp.l;
            type = hp.type;
            s = hp.scale;
            switch type
                case 'squared exponential iso'
                    if length(L) == 1
                        obj.kernel = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(L^2)));
                    else
                        error('Lengthscale parameter should be scalar!');
                    end
                case 'linear iso'
                    if isempty(L)
                        obj.kernel = @(x1,x2) s * x1(:)'*x2(:);
                    else
                        error('Lengthscale parameter should be empty!');
                    end
                case 'squared exponential ard'
                    InvGamma = diag(1./(L.^2));
                    obj.kernel = @(x1,x2) s * exp(-0.5*((x1(:)-x2(:))')*InvGamma*(x1(:)-x2(:)));
                case 'linear ard'
                    InvGamma = diag(1./(L.^2));
                    obj.kernel = @(x1,x2) s*x1(:)'*InvGamma*x2(:);
                otherwise
                    error('Unrecognized kernel type');
            end
        
        end        
        
        % Constructs the kernel matrix all at once.
        % Only used to generate test functions.
        % INPUTS: 
        % x - training data points
        % kernel - covariance function
        % OUTPUTS
        % out - matrix of size nxn, where n is the length of us
        function buildCov(obj,x)

            len = size(x,2);
            Kmat = zeros(len,len);
            for i = 1:len
                for j = i:len
                    Kmat(i,j) = obj.kernel(x(:,i), x(:,j));
                end
            end
            obj.cov = Kmat + Kmat' - diag(diag(Kmat));
            obj.cov = obj.cov + obj.hp.noise.var * eye(len);
        end
        
        % Construct the kernel matrix iteratively
        % INPUTS: 
        % t - last time stage (up to horizon)
        % us - control law points 
        % fun - function handle for the kernel
        % mat - previous kernel matrix
        % OUTPUT
        % Kmat - updated kernel matrix at the previous points and current test pt.
        function updateCov(obj,xstar)

            % updates iteratively
            Kmat = obj.cov;
            k = obj.covTestAndData(xstar);
            % since kernel matrix is symmetric 
            % simply update the last row with 
            % the last column
            Kmat(:,end+1) = k;
            kxx = obj.kernel(xstar,xstar) + obj.hp.noise.var;
            Kmat(end+1,:) = [k;kxx]';
            obj.cov = Kmat;
        end
        
        % construct the kernel vector - covariances of the test pt and previous pts
        % INPUTS: 
        % xstar - test point
        % x     - data points 
        % OUTPUT
        % cov - kernel column vector describing the covariances between
        %       the previous points and current test pt.
        %
        function cov = covTestAndData(obj,xstar)

            n = size(obj.x,2);
            cov = zeros(n,1);
            for i = 1:n
                cov(i) = obj.kernel(obj.x(:,i),xstar);
            end
        
        end
        
        
    end

end