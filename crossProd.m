% Crossproduct
% Gives the unit vector that is orthogonal to both unit vectors
% Inputs must be orthogonal

function out = crossProd(v1,v2)

    eps = 1e-3;
    v1 = v1(:);
    v2 = v2(:);

    assert(length(v1) == 3);
    assert(length(v2) == 3);
    % make sure v1 and v2 are orthogonal
    assert(v1'*v2 < eps,'Inputs must be orthogonal');
    
    % make sure v1 and v2 are unit norm    
    if v1'*v1 ~= 1.000
        % give warning and normalize
        warning('First input not normalized! Normalizing...');
        v1 = v1./norm(v1,2);
    end
    
    if v2'*v2 ~= 1.000
        warning('Second input not normalized! Normalizing...');
        v2 = v2./norm(v2,2);
    end
    
    out(1) = v1(2)*v2(3) - v1(3)*v2(2);
    out(2) = v1(3)*v2(1) - v1(1)*v2(3);
    out(3) = v1(1)*v2(2) - v1(2)*v2(1);
    
    % make sure out is normalized
    out = out(:);
    out = out./norm(out,2);
    
    assert(out'*v1 < eps);
    assert(out'*v2 < eps);
    


end