function [R1,N1]=downsample_index(n,b)

num_t = 2 ^ b;
num_f = 2 ^ (n-b);

%index in the range space of downsample matrix
for i = 1 : num_t
    R1(i) = (i-1); 
end
%index in the null space of downsample matrix
for i = 1 : num_f
    N1(i) = (i-1) * num_t;
end

end


