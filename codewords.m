%all possible code word for different cases.
function [CW1,CW2,CW3] = codewords(n_alias)
%n_alias = 8;

%choose the shift in the order
H = hadamard(n_alias);
order = order_shift(n_alias);
H_reorder = zeros(n_alias,n_alias);
for i = 1:n_alias
    H_reorder(i,:)=H(order(i),:);
end

%all possible cases
C1 = cases_order1(n_alias);
C2 = cases_order2(n_alias);
C3 = cases_order3(n_alias);

%generate the codewords
CW1 = C1*H_reorder';
CW2 = C2*H_reorder';
CW3 = C3*H_reorder';
end


