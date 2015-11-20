bottom_layer_size=100; % number of nodes on bottom layer
filename=strcat(int2str(bottom_layer_size),'pyramid.txt');
fileID=fopen(filename,'w');

totalnodes = bottom_layer_size*(bottom_layer_size+1)/2;
totalelements = bottom_layer_size*(bottom_layer_size-1);



%prints number of elements and nodes
fprintf(fileID,'%d, %d\n',totalelements,totalnodes);


%prints elements (node1,node2)
start=1;
current_size=bottom_layer_size;

for j=1:bottom_layer_size-1
    ending=start+current_size-1;
    makelayer(fileID,start,ending);
    start=start+current_size;
    current_size = current_size-1;

end

%prints constraint for bottom layer
for i=1:bottom_layer_size
    fprintf(fileID,'%d, %d\n%d, %d\n%d, %d\n',i,1,i,2,i,3);
end
%prints constraint for the rest (z direction fixed)
for i=bottom_layer_size+1:totalnodes
    fprintf(fileID,'%d, %d\n',i,3);
end

%constraint done flag
fprintf(fileID,'0, 0\n');
xoffset=10;
yoffset=5;
for j=0:bottom_layer_size-1
    for i=0:bottom_layer_size-1-j
        xcoor = i*xoffset+j*xoffset/2;
        ycoor=j*yoffset;
        zcoor=0;
        fprintf(fileID,'%d, %d, %d\n',xcoor,ycoor,zcoor);
    end
end

%material property
for i=1:totalelements
   area=1;
   modulus=29000;
   fprintf(fileID,'%d, %d\n',area,modulus);
end


for i=bottom_layer_size+1:totalnodes
    direction=2;
    load=-10;
    fprintf(fileID,'%d, %d, %d\n',i,direction,load);
end

fprintf(fileID,'0, 0, 0\n');
fprintf(fileID,'2, 0.1, 0.1 \n200 \n0.001, 0.001');

fclose(fileID);