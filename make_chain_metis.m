elements = 999; %number of elements: elements should be odd for a truss bridge
nodes=(elements+1)/2+1;  % number of nodes
weight= 1;
filename = strcat(int2str(elements),'elementschain.txt');
fileID=fopen(filename,'w');


fprintf(fileID,'%d, %d\n', elements, weight); %% write #element, weight

%the diagonal truss members
diagonal= (elements+1)/2;

for i=1:diagonal
    fprintf(fileID,'%d, %d\n', i, i+1); %% node 1&2 of element
end

%truss members at the bottom
bottom= (elements-diagonal+1)/2;
for i=1:2:bottom
   fprintf(fileID, '%d, %d\n', i, i+2);
end

%truss members at the top
top= bottom-1;
for i=1:2:top
   fprintf(fileID, '%d, %d\n', i+1, i+3);
end

fclose(fileID);
