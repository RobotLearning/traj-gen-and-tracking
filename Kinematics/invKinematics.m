%% General Inverse Kinematics Function to be used for BarrettWAM

% Ahmats are the homogenous matrices of the last joint 
% 
% Taken and adapted from The Robotics Toolbox for Matlab (RTB) by
% Guilherme Maeda, 2015.
% Copyright (C) 1993-2011, by Peter I. Corke
% http://www.petercorke.com

function qTotal = invKinematics(Ahmats, varargin)

    %  set default parameters for solution
    opt.ilimit = 1000;
    opt.tol = 1e-6;
    opt.alpha = 1;
    opt.plot = false;
    opt.pinv = false;
    opt.varstep = true;

    [opt,args] = tb_optparse(opt, varargin);
    q0 = args{1};
    q = q0;
    numDofs = size(q0,1); % num of dofs

    % robot.ikine(tr, q, m)
    m = args{2};
    m = m(:);
    if numel(m) ~= 6
        error('Mask matrix should have 6 elements');
    end
    if numel(find(m)) ~= numDofs
        error('Mask matrix must have same number of 1s as robot DOF')
    end

    % make this a logical array so we can index with it
    m = logical(m);

    numPoints = size(Ahmats,3);    % number of points
    qTotal = zeros(numPoints, numDofs);  % preallocate space for results
    tcount = 0; % total iteration count

    if ~ishomog(Ahmats)
        error('RTB:ikine:badarg', 'T is not a homog xform');
    end

    % compute jacobian of q0 to check for initial conditioning
    [xLink,xOrigin,xAxis,~] = barrettWamKinematics(q0,obj.PAR);
    J0 = jacobian(xLink,xOrigin,xAxis);

    if cond(J0) > 100
        warning('RTB:ikine:singular', ...
        'Initial joint angles results in near-singular configuration, this may slow convergence');
    end

    for i = 1:numPoints
        
        T = Ahmats(:,:,i);

        normErr = Inf;
        % initialize state for the inverse kinematics loop
        % previous error
        eprev = Inf;
        save.e = [Inf Inf Inf Inf Inf Inf];
        save.q = [];
        count = 0;

        while true
            % update the count and test against iteration limit
            count = count + 1;
            if count > opt.ilimit
                warning('ikine: iteration limit %d exceeded (row %d), final err %f', ...
                    opt.ilimit, i, normErr);
                q = NaN*ones(1,numDofs);
                break
            end

            % compute the error
            [xLink,xOrigin,xAxis,~] = barrettWamKinematics(q,obj.PAR);
            x = xLink(6,:)';
            e = tr2delta(x, T);

            % optionally adjust the step size
            if opt.varstep
                % test against last best error, only consider the DOF of
                % interest
                if norm(e(m)) < norm(save.e(m))
                    % error reduced,
                    % let's save current state of solution and rack up the step size
                    save.q = q;
                    save.e = e;
                    opt.alpha = opt.alpha * (2.0^(1.0/8));
                    if opt.verbose > 1
                        fprintf('raise alpha to %f\n', opt.alpha);
                    end
                else
                    % rats!  error got worse,
                    % restore to last good solution and reduce step size
                    q = save.q;
                    e = save.e;
                    opt.alpha = opt.alpha * 0.5;
                    if opt.verbose > 1
                        fprintf('drop alpha to %f\n', opt.alpha);
                    end
                end
            end

            % compute the Jacobian
            J = jacobian(xLink,xOrigin,xAxis);
            %J = jacob0(obj, q);

            % compute change in joint angles to reduce the error, 
            % based on the square sub-Jacobian
            if opt.pinv
                dq = opt.alpha * pinv(J(m,:)) * e(m);
            else
                dq = J(m,:)' * e(m);
                dq = opt.alpha * dq;
            end

            % diagnostic stuff
            if opt.verbose > 1
                fprintf('%d:%d: |e| = %f\n', i, count, normErr);
                fprintf('       e  = '); disp(e');
                fprintf('       dq = '); disp(dq');
            end

            % update the estimated solution
            q = q + dq';
            normErr = norm(e(m));

            if norm(e) > 1.5 * norm(eprev)
                warning('RTB:ikine:diverged', 'solution diverging, try reducing alpha');
            end
            eprev = e;

            if normErr <= opt.tol
                break
            end

        end  % end ikine solution for tr(:,:,i)
        
        qTotal(:,i) = q';
        tcount = tcount + count;
        if opt.verbose
            fprintf('%d iterations\n', count);
        end
    end
    
    if opt.verbose && numPoints > 1
        fprintf('TOTAL %d iterations\n', tcount);
    end
    
end