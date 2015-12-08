function [data_new, data_test ]= remove_test_data(data, testSet)

    nJoints = size(data,2);
    index = 1;
    
    for j = 1:nJoints
        data_new(j).q       = [];
        data_new(j).qdot    = [];
    end
    
    
    for d = 1:size(data(1).q, 1)
        
        for j=1:nJoints

            if d == testSet
                data_test(j).q      = data(j).q(d,:);
                data_test(j).qdot   = data(j).qdot(d,:);
            else
                data_new(j).q    = [data_new(j).q;    data(j).q(d,:)];
                data_new(j).qdot = [data_new(j).qdot; data(j).qdot(d,:)];
            end
        
        end        
    end
    

end

