% Calculate offset between camera 1 and camera 3

function offset = calculateOffset(t1,b1,t3,b3)

idx1And3 = zeros(length(t1),1);
idx3And1 = zeros(length(t3),1);
for i = 1:length(t1)
    for j = 1:length(t3)
        if abs(t1(i) - t3(j)) < 1e-3
            idx1And3(i) = 1;
            idx3And1(j) = 1;
        end
    end
end

N = sum(idx1And3);
idx1And3 = logical(idx1And3);
idx3And1 = logical(idx3And1);

b1 = b1(idx1And3,:);
b3 = b3(idx3And1,:);

offset = (1/N) * sum(b1 - b3,1);

end