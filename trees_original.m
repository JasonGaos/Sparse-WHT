%generate the tree structure for decoder
function [tree_1,tree_2,tree_3] = trees(CW1,CW2,CW3,num_co)
%num_co = 8;
[CW1,CW2,CW3]=codewords(num_co);
%tree for order 1 
[r,c]=size(CW1);
node_max = 2*num_co;
tree_1 = cell(node_max,1);
for i =1:r
    p=1;
    for j=1:log2(c)+1
        tree_1{p,1} = [tree_1{p,1},i];
        if CW1(i,j+1) == 1
            p = 2 * p;
        elseif CW1(i,j+1) == -1
            p = 2 * p + 1;
        end               
    end
end
% L1 = zeros(r,1);
% for i = 1:node_max
%     if length(tree_1{i}) == 1 && L1(tree_1{i}) == 0
%         L1(tree_1{i})=length(dec2base(i,2));
%     end
% end

%tree for order 2 
[r,c]=size(CW2);
node_max = (3^num_co-1)/2;
i = 0;
while 3^i < node_max
    level_2(i+1)=3^i;
    i = i+1;
end
level_2 = cumsum(level_2);
tree_2 = cell(node_max,1);
for i =1:r
    p=1;
    tree_2{p,1} = [tree_2{p,1},i];
    for j=1:c-1
        if CW2(i,j+1) == 2
            p = 3 * p - 1;
        elseif CW2(i,j+1) == 0
            p = 3 * p;
        elseif CW2(i,j+1) == -2
            p = 3 * p + 1;
        end               
        tree_2{p,1} = [tree_2{p,1},i];
    end
end
% L2 = zeros(r,1);
% for i = 1:node_maxv
%     if length(tree_2{i}) == 1
%         for j = 1:length(level_2)
%             if i<=level_2(j) && L2(tree_2{i}) == 0
%                 L2(tree_2{i})=j;
%             end
%         end
%     end
% end

%tree for order 3 

[r,c]=size(CW3);
node_max = (4^num_co-1)/3;
i = 0;
while 4^i < node_max
    level_3(i+1)=4^i;
    i = i+1;
end
level_3 = cumsum(level_3);
tree_3 = cell(node_max,1);
for i =1:r
    p=1;
    tree_3{p,1} = [tree_3{p,1},i];
    for j=1:c-1
        if CW3(i,j+1) == 3
            p = 4 * p - 2;
        elseif CW3(i,j+1) == 1
            p = 4 * p -1;
        elseif CW3(i,j+1) == -1
            p = 4 * p;
        elseif CW3(i,j+1) == -3
            p = 4 * p + 1;
        end
        tree_3{p,1} = [tree_3{p,1},i];
    end
end
% L3 = zeros(r,1);
% for i = 1:node_max
%     if length(tree_3{i}) == 1
%         for j = 1:length(level_3)
%             if i<=level_3(j) && L3(tree_3{i}) == 0
%                 L3(tree_3{i})=j;
%             end
%         end
%     end
% end

% avg_length(1) = mean(L1);
% avg_length(2) = mean(L2);
% avg_length(3) = mean(L3);

end

