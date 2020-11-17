function [R1,N1,R2,N2] = downsample_index_nonoverlap(n,b)

num_t = 2 ^ b; 
num_f = 2 ^ (n-b);

%index in the range space of downsample matrix 1
for i = 1 : num_t
    R1(i) = (i-1); 
end
%index in the null space of downsample matrix 1
for i = 1 : num_f
    N1(i) = (i-1) * num_t;
end


%index in the range space of downsample matrix 2
for i = 1 : num_t
    R2(i) = (i-1) * num_f; 
end
%index in the null space of downsample matrix 2
for i = 1 : num_f
    N2(i) = (i-1);
end



end