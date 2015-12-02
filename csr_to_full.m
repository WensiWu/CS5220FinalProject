function [ full] = csr_to_full(a, ia , ja )
n = length(ia)-1;
full = zeros(n);
size = ia(n+1)-1;

for i=1:n
    number_in_row=ia(i+1)-ia(i);
    for j=0:number_in_row-1
        full(i,ja(ia(i)+j))=a(ia(i)+j);
    end
end

        

end

