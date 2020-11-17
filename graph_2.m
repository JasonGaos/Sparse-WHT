%%
%parameters
performance = zeros(10,1);
for sparsity = 1:10
n_v = 8;
n_e = sparsity;
b1 = n_v-3;
b2 = n_v-3;
N = 2^n_v;
L1 = order_shift(2 ^ (n_v - b1));
L2 = order_shift(2 ^ (n_v - b2));
n=100;
extra = 0;
success = 0;

for test_round = 1:1000
%% generate data
W = zeros(1,2^n_v);


index_possible = graph_index(n_v);
index_choosen = [];

while length(index_choosen)<sparsity
a = rand;
i = round(a*((n_v*(n_v-1)/2 -1))) + 1;
index_choosen = [index_choosen;index_possible(i)];
index_choosen = intersect(index_choosen,index_possible);
end


for i = 1:sparsity
    W(index_choosen(i)) = 1;
end


%W = [ 1,0,0,0, 0,0,1,0, 0,0,0,0 ,0,0,0,0];
w = ifwht(W,length(W),'hadamard');
W1 = W;

%% down-sample
% [R1,N1]=downsample_index(n_v,b1);
% [R2,N2]=downsample_index(n_v,b2);
[R1,N1,R2,N2] = downsample_index_nonoverlap(n_v,b1);
%tree for stage 1
[CW11,CW12,CW13] = codewords(2^(n_v-b1));
[tree_11,tree_12,tree_13] = trees(CW11,CW12,CW13,2^(n_v-b1),extra);
%all possible cases
C11 = cases_order1(2^(n_v - b1));
C12 = cases_order2(2^(n_v - b1));
C13 = cases_order3(2^(n_v - b1));
% order of shift
order1 = order_shift(2^(n_v - b1));

%tree for stage 2
[CW21,CW22,CW23] = codewords(2^(n_v-b2));
[tree_21,tree_22,tree_23] = trees(CW21,CW22,CW23,2^(n_v - b2),extra);
%all possible cases
C21 = cases_order1(2^(n_v - b2));
C22 = cases_order2(2^(n_v - b2));
C23 = cases_order3(2^(n_v - b2));
% order of shift
order2 = order_shift(2^(n_v-b2));

stage1_indicator = 1;
stage2_indicator = 1;
round_peeling = 0;

while (W*W'~=0)&& round_peeling<=n_e 
    round_peeling = round_peeling +1;
%% STAGE 1
X_recover = zeros(1,2^n_v);
 
number_shift = 0; %starts from one
shift_indicator = 1;
stage2_indicator = 0;
type_tree = zeros(2^(n_v - b1),1);
pointer = ones(length(R1),1);%point the index of the tree, starts from one
index = [];
shift = 0;

%% peeling decoder

while shift_indicator == 1 && shift<=n_v-b1%while
    %new shift downsample and shorter WHT
    number_shift = number_shift + 1;
    shift = shift + 1;
    if shift >= (n_v-b1)+2+extra
        break
    end
    x1 = [];
    for i = 1:length(R1)
        x1=[x1,w(R1(i)+1+N1(order1(shift)))];
    end
    X1 = fwht(x1,length(x1),'hadamard');
for i = 1:length(X1)
    %initial the type_tree
    if number_shift == 1
        type_tree(i) = X1(i);
    end
    %pointer changes in type 1    
    if type_tree(i) == 1
        if number_shift ~= 1
            if X1(i) == 1
                pointer(i) = pointer(i) * 2;
            elseif X1(i) == -1
                pointer(i) = pointer(i) * 2 + 1;
            end
        end        
        a = size(tree_11{pointer(i)});
        if a == [1,1]
            temp = tree_11(pointer(i));
            distribution  = C11(temp{1},:);
            k = find(distribution);
            index = [index,(k-1)*2^b1 + i];
            shift_indicator = 0;
        else
            shift_indicator = 1;
        end
        %pointer changes in type 2            
    elseif type_tree(i) == 2
        if number_shift ~= 1
            if X1(i) == 2
                pointer(i) = 3 * pointer(i) - 1;
            elseif X1(i) == 0
                pointer(i) = 3 * pointer(i);
            elseif X1(i) == -2
                pointer(i) = 3 * pointer(i) + 1;
            end
        end
        a = size(tree_12{pointer(i)});
        if a == [1,1]
            temp = tree_12(pointer(i));
            distribution  = C12(temp{1},:);
            k = find(distribution);
            index = [index,(k-1)*2^b1 + i];
            shift_indicator = 0;
        else
            shift_indicator = 1;
        end
    end
end
end
if size(index)~=[0,0]
    for i = 1 : length(index)
        X_recover(index(i))=1;
    end
end

w = w - ifwht(X_recover,length(X_recover),'hadamard');
%% STAGE 2
X_recover = zeros(1,2^n_v);

number_shift = 0; %starts from one
shift_indicator = 1;
stage1_indicator = 1;
type_tree = zeros(2^(n_v - b2),1);
pointer = ones(length(R2),1);%point the index of the tree, starts from one
index = [];
shift = 0;

%% peeling decoder

while shift_indicator == 1 && shift<=n_v-b2%while
    %new shift downsample and shorter WHT
    number_shift = number_shift + 1;
    shift = shift + 1;
    if shift >= n_v-b2+2+extra
        break
    end
    x2 = [];
    for i = 1:length(R2)
        x2=[x2,w(R2(i)+1+N2(order2(shift)))];
    end
    X2 = fwht(x2,length(x2),'hadamard');
for i = 1:length(X2)
    %initial the type_tree
    if number_shift == 1
        type_tree(i) = X2(i);
    end
    %pointer changes in type 1    
    if type_tree(i) == 1
        if number_shift ~= 1
            if X2(i) == 1
                pointer(i) = pointer(i) * 2;
            elseif X2(i) == -1
                pointer(i) = pointer(i) * 2 + 1;
            end
        end        
        a = size(tree_21{pointer(i)});
        if a == [1,1]
            temp = tree_21(pointer(i));
            distribution  = C21(temp{1},:);
            k = find(distribution);
            index = [index,(k-1)*2^b2 + i];
            shift_indicator = 0;
        else
            shift_indicator = 1;
        end
        %pointer changes in type 2            
    elseif type_tree(i) == 2
        if number_shift ~= 1
            if X2(i) == 2
                pointer(i) = 3 * pointer(i) - 1;
            elseif X2(i) == 0
                pointer(i) = 3 * pointer(i);
            elseif X2(i) == -2
                pointer(i) = 3 * pointer(i) + 1;
            end
        end
        a = size(tree_22{pointer(i)});
        if a == [1,1]
            temp = tree_22(pointer(i));
            distribution  = C22(temp{1},:);
            k = find(distribution);
            index = [index,(k-1)*2^b2 + i];
            shift_indicator = 0;
        else
            shift_indicator = 1;
        end
    end
end
end

if size(index)~=[0,0]
    for i = 1 : length(index)
        X_recover(index(i))=1;
    end
end

w = w - ifwht(X_recover,length(X_recover),'hadamard');


%%check
W = fwht(w,length(w),'hadamard');
end

%%succuss rate

if W*W' == 0
success = success+1;
end
performance(sparsity)=success;
end
end