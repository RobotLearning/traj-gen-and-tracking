% Checks if the point is inside a 2d polygon

function flag = checkIfPointIsInsidePolygon(polygon,p)
  
    counter = 0;

    N = size(polygon,2);
    p1 = polygon(1,:);
    X = 1; % indices
    Y = 2;
    
    for i = 2:N
        p2 = polygon(i,:);
        if p(Y) > min(p1(Y),p2(Y)) 
            if p(Y) <= max(p1(Y),p2(Y)) 
                if p(X) <= max(p1(X),p2(X)) 
                    if p1(Y) ~= p2(Y) 
                        xinters = (p(Y)-p1(Y))*(p2(X)-p1(X))/(p2(Y)-p1(Y))+p1(X);
                        if p1(X) == p2(X) || p(X) <= xinters
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
        p1 = p2;
    end
        
    flag = rem(counter,2);
end
