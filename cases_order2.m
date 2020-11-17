%all possible cases contains two non-zero coefficients
function A = cases_order2(n)
d = 2;
number = nchoosek(n,d);
A = zeros(number,n);
num = 1;
while num <= number
    for i = 1:n-1
        for j = i+1:n
                A(num,i) = 1;
                A(num,j) = 1;
                num = num + 1;
        end
    end
end

end