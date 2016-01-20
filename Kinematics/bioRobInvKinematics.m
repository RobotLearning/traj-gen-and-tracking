%% Inverse Kinematics Function for the BioRob

% TODO: replace the jacobian and implement forward kinematics

% Taken and adapted from The Robotics Toolbox for Matlab (RTB) by
% Guilherme Maeda, 2015.
% Copyright (C) 1993-2011, by Peter I. Corke
% http://www.petercorke.com

function [qt, histout] = bioRobInvKinematics(tr, varargin)
    %  set default parameters for solution
    opt.ilimit = 1000;
    opt.tol = 1e-6;
    opt.alpha = 1;
    opt.plot = false;
    opt.pinv = false;
    opt.varstep = true;

    [opt,args] = tb_optparse(opt, varargin);

    numDofs = 4; % num of dofs

    % robot.ikine(tr, q)
    if ~isempty(args)
        q = args{1};
        if isempty(q)
            q = zeros(1, numDofs);
        else
            q = q(:)';
        end
    else
        q = zeros(1, numDofs);
    end

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

    npoints = size(tr,3);    % number of points
    qt = zeros(npoints, numDofs);  % preallocate space for results
    tcount = 0;              % total iteration count

    if ~ishomog(tr)
        error('RTB:ikine:badarg', 'T is not a homog xform');
    end

    J0 = jacob0(robot, q);
    
    % J0 = J0(m, m); %<================= WHY DO I HAVE TO DO THIS????

    if cond(J0) > 100
        warning('RTB:ikine:singular', 'Initial joint angles results in near-singular configuration, this may slow convergence');
    end

    history = [];
    failed = false;
    for i = 1:npoints
        T = tr(:,:,i);

        nm = Inf;
        % initialize state for the ikine loop
        eprev = Inf;
        save.e = [Inf Inf Inf Inf Inf Inf];
        save.q = [];
        count = 0;

        while true
            % update the count and test against iteration limit
            count = count + 1;
            if count > opt.ilimit
                warning('ikine: iteration limit %d exceeded (row %d), final err %f', ...
                    opt.ilimit, i, nm);
                failed = true;
                q = NaN*ones(1,numDofs);
                break
            end

            % compute the error
            e = tr2delta(robot.fkine(q'), T);

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
            J = jacob0(robot, q);

            % compute change in joint angles to reduce the error, 
            % based on the square sub-Jacobian
            if opt.pinv
                dq = opt.alpha * pinv( J(m,:) ) * e(m);
            else
                dq = J(m,:)' * e(m);
                dq = opt.alpha * dq;
            end

            % diagnostic stuff
            if opt.verbose > 1
                fprintf('%d:%d: |e| = %f\n', i, count, nm);
                fprintf('       e  = '); disp(e');
                fprintf('       dq = '); disp(dq');
            end

            % update the estimated solution
            q = q + dq';
            nm = norm(e(m));

            if norm(e) > 1.5*norm(eprev)
                warning('RTB:ikine:diverged', 'solution diverging, try reducing alpha');
            end
            eprev = e;

            if nm <= opt.tol
                break
            end

        end  % end ikine solution for tr(:,:,i)
        qt(i,:) = q';
        tcount = tcount + count;
        if opt.verbose
            fprintf('%d iterations\n', count);
        end
    end
    
    if opt.verbose && npoints > 1
        fprintf('TOTAL %d iterations\n', tcount);
    end
    
end