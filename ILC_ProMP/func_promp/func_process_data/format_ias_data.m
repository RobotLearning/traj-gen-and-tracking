function data_ias = format_ias_data(data)
% This is simply to put the data on the format that is used by the IAS
% promp class

    [nJoints] = size(data,2);
    [nDemo] = size(data(1).q);
    
    for k=1:nDemo
        
        collect = [];
        
        for j=1:nJoints
            collect = [collect data(j).q(k,:)'];
        end
        data_ias{k,:} = collect;
        %data_ias{k,:} = [data(1).q(k,:)' data(2).q(k,:)'...
        %                 data(3).q(k,:)' data(4).q(k,:)'];
    end

end