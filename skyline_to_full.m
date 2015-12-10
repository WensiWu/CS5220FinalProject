function [ full ] = skyline_to_full(a,max)

n = length(max)-1;
full = sparse(length(max)-1, length(max)-1);
    
for j=1:n
    number_in_col = max(j+1)- max(j);
    for i = 0:number_in_col-1
        full(j-number_in_col+i+1,j)=a(max(j)+i);
    end
end

end

