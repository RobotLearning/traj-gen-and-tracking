%% General Inverse Kinematics Function to be used for BarrettWAM

% Ahmats are the homogenous matrix of the desired end effector pose 
% 
% Taken and adapted from The Robotics Toolbox for Matlab (RTB) by
% Guilherme Maeda, 2015
% Okan Koc, 2016.
%
% Copyright (C) 1993-2011, by Peter I. Corke
% http://www.petercorke.com
%
function qTotal = invKinematics(Ahmats, q0, PAR)

    %  set default parameters for solution
    opt.ilimit = 1000;
    opt.tol = 1e-6;
    opt.alpha = 1;
    opt.plot = false;
    opt.pinv = false;
    opt.varstep = true;
    opt.verbose = false;

    q = q0;
    numDofs = size(q0,1); % num of dofs
    m = ones(1,6); % not masking any variables as of now 

    % make this a logical array so we can index with it
    m = logical(m);

    numPoints = size(Ahmats,3);    % number of points
    qTotal = zeros(numDofs, numPoints);  % preallocate space for results
    tcount = 0; % total iteration count

    if ~ishomog(Ahmats)
        error('T is not a homog xform');
    end

    % compute jacobian of q0 to check for initial conditioning
    [xLink,xOrigin,xAxis,~] = barrettWamKinematics(q0,PAR);
    J0 = jacobian(xLink,xOrigin,xAxis);

    if cond(J0) > 100
        warning('Initial joint angles results in near-singular configuration, this may slow convergence');
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
                warning('InvKin: iteration limit %d exceeded (row %d), final err %f', ...
                    opt.ilimit, i, normErr);
                %q = NaN*ones(1,numDofs);
                break;
            end

            % compute the error
            [xLink,xOrigin,xAxis,As] = barrettWamKinematics(q,PAR);
            Tnow = squeeze(As(6,:,:));
            e = getDiffBetweenPoses(Tnow, T);

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
            q = q + dq;
            normErr = norm(e(m));

            if norm(e) > 1.5 * norm(eprev)
                warning('solution diverging, try reducing alpha');
            end
            eprev = e;

            if normErr <= opt.tol
                break
            end

        end  % end ikine solution for tr(:,:,i)
        
        qTotal(:,i) = q;
        tcount = tcount + count;
        if opt.verbose
            fprintf('%d iterations\n', count);
        end
    end
    
    if opt.verbose && numPoints > 1
        fprintf('TOTAL %d iterations\n', tcount);
    end
    
end

% checks to see if the homogenous matrix tr is valid
function h = ishomog(tr)

    d = size(tr);
    eps = 1e-3;
    
    if ndims(tr) >= 2
        h =  all(d(1:2) == [4 4]);

        if h 
            h = abs(det(tr(1:3,1:3)) - 1) < eps;
        end

    else
        h = false;
    end
    
end

% get differential difference between poses T0 and T1
function delta = getDiffBetweenPoses(T0, T1)

    eps = 1e-3;
    R0 = T0(1:3,1:3); 
    assert(abs(det(R0) - 1) < eps,'T0 is not valid!');
    R1 = T1(1:3,1:3);
    assert(abs(det(R1) - 1) < eps,'T1 is not valid!');
    % in world frame
    vec = T1(1:3,4) - T0(1:3,4);
    delta = [vec; skewSymToVec(R1*R0'-eye(3,3))];

end

% convert skew symmetric matrix to vector
% no checking done
function vec = skewSymToVec(S)

vec = 0.5*[S(3,2)-S(2,3); 
           S(1,3)-S(3,1); 
           S(2,1)-S(1,2)];

end