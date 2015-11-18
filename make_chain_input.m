elements = 2; % number of elements
nodes=elements+1; % number of nodes
filename = strcat(int2str(elements),'elementschain.txt');
fileID=fopen(filename,'w');


fprintf(fileID,'%d, %d\n', elements, nodes); %% write #element, #node

for i=1:elements
    fprintf(fileID,'%d, %d\n', i, i+1); %% node 1&2 of element
end

%joint constraints
for i=1:nodes 
    if mod(i,2)==1
        for j=1:3
            fprintf(fileID,'%d, %d\n', i, j);
        end
    else 
        fprintf(fileID,'%d, %d\n', i, 3);
    end
end

%done with constraint statement:
fprintf(fileID,'0, 0\n');


xoffset=10;
yoffset=2;
for i=0:nodes-1
    xcoor=i*xoffset;
    ycoor=mod(i,2)*yoffset;
    zcoor=0;
    fprintf(fileID,'%d, %d, %d\n', xcoor, ycoor, zcoor); %% coordinate of node i
end

% assigns material props to elements
for i=1:elements
    area = 1; 
    modulus = 29000;
    fprintf(fileID,'%d, %d\n', area, modulus); %% assign material props : area and modulus
end

% assigns load to even numbered node
for i=1:nodes
    if mod(i,2)==0
         load = 10;
         fprintf(fileID,'%d, %d, %d\n', i, 2, -load);
    end
end
%signals done with loading   
fprintf(fileID,'0, 0, 0\n');

fclose(fileID);
