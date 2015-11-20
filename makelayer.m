function [  ] = makelayer(fileID,a,b)

% n = length of current layer
n = b-a+1;

fprintf(fileID,'%d, %d\n',a,a+n);
for i=1:n-2
    fprintf(fileID,'%d, %d\n',a+i, a+i+n-1);
    fprintf(fileID,'%d, %d\n',a+i, a+i+n);
end
fprintf(fileID,'%d, %d\n',b,b+n-1);

end

