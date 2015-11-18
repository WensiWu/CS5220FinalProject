elements = 999; % number of elements
nodes=elements+1; % number of nodes
filename = strcat(int2str(elements),'elementschain.txt');
fileID=fopen(filename,'w');


fprintf(fileID,'%d, %d\n', elements, nodes); %% write #element, #node

for i=1:elements
    fprintf(fileID,'%d, %d\n', i, i+1); %% node 1&2 of element
end

fprintf(fileID,'0, 0\n');

xoffset=10;
yoffset=2;
for i=1:nodes
    xcoor=i*xoffset;
    ycoor=rem(i+1,2)*yoffset;
    zcoor=0;
    fprintf(fileID,'%f, %f, %f\n', xcoor, ycoor, zcoor); %% coordinate of node i
end

% assigns material props to elements
for i=1:elements
    area = randi([1,4]); 
    modulus = 20000+10000*rand;
    fprintf(fileID,'%f, %f\n', area, modulus); %% assign material props : area and modulus
end

% assigns load to even numbered node
for i=1:nodes
    if rem(i,2)==0
         load = 5+20*rand;
         fprintf(fileID,'%d, %d, %f\n', i, 2, -load);
    end
end
%signals done with loading   
fprintf(fileID,'0, 0, 0\n');

fclose(fileID);
