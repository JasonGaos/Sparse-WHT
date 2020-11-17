function [A,B]=downsample_matrix(n,downsample_ratio)
% downsample_ratio = 0.3;
% n=8;
b = ceil(downsample_ratio * n);

downsample_matrix_1 = [eye(b),zeros(b,n - b)];
downsample_matrix_2 = [zeros(b,n - b),eye(b)];

A = downsample_matrix_1; 
B = downsample_matrix_2;
end