% Checks if the point is inside a 2d polygon

function cross = checkIfPointIsInsidePolygon(polygon,p)
  
    N = size(polygon,1);    
    cross = false;

    f = zeros(N,1);
    beq = [p(:);1];
    A = polygon';
    Aeq = [A;ones(1,N)];
    lb = zeros(N,1);
    ub = ones(N,1);
    [x,fval,exitflag,output] = linprog(f,[],[],Aeq,beq,lb,ub);
    
    if exitflag ~= -2
        % solution was found
        cross = true;
    end

%     p1 = polygon(1,:);
%     X = 1; % indices
%     Y = 2;
%     for i = 2:N
%         p2 = polygon(i,:);
%         if (p2(Y) > p(Y)) ~= (p1(Y) > p(Y)) && ...
%             (p(X) < (p1(X)-p2(X))*(p(Y)-p2(Y))/(p1(Y)-p2(Y)) + p2(X))
%             cross = ~cross;
%         end
%         p1 = p2;
%     end
%     

end
