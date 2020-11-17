%order of shifts. gives the index of the row from a hadamard matrix of order n_v
function L = order_shift(n_v)
L = zeros(n_v,1);
L(1) = 1;
count = 2;
j = n_v+1;
while j>2
    j = (j+1)/2;
    L(count) = j;
    count = count+1;
end

for i = 1:n_v
    p = ismember(i,L);
    if p == 0
        L(count) = i;
        count = count + 1;
    end
end
end


