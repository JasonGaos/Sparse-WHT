%all possible cases contains three non-zero coefficients
function A = cases_order3(n)
d = 3;
number = nchoosek(n,d);
A = zeros(number,n);
num = 1;
while num <= number
    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n 
                A(num,i) = 1;
                A(num,j) = 1;
                A(num,k) = 1;
                num = num + 1;
            end
        end
    end
end

end