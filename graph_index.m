%generate a set of all possible index for graph data
%the binary representation contains two '1'.

function b = graph_index(n)
% n = 8;
index_binary = cases_order2(n);
num = n*(n-1)/2;
b = zeros(num,1);

for i = 1:num
    for j = 1:n
        b(i) = b(i) + index_binary(i,j)*2^(j-1);
    end
end

end