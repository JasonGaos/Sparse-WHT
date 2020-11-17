function [R1,N1,R2,N2]=downsample_index2(n,b1,b2)
% n = 4;
% b1 = 1;
% b2 = 2;
% 
num_t1 = 2^b1;
num_t2 = 2^b2;
num_f1 = 2^(n - b1);
num_f2 = 2^(n - b2);

%index in the range space of downsample stage 1
for i = 1 : num_t1
    R1(i) = i-1;
end

for i = 1: num_t1
    for j = 1 : num_f1
        N1(i,j) = (i-1)*num_f1 + j;
    end
end

%index in the range space of downsample stage 2
for i = 1 : num_t2
    R2(i) = i-1;
end

for i = 1: num_t2
    for j = 1 : num_f2
        N2(i,j) = (i-1)*num_f2 + j;
    end
end

end