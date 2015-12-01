elements = 99999; %number of elements: elements should be odd for a truss bridge
nodes=(elements+1)/2+1;  % number of nodes
filename = strcat(int2str(elements),'elementschain.txt');
fileID=fopen(filename,'w');


fprintf(fileID,'%d, %d\n', elements, nodes); %% write #element, #node

%the diagonal truss members
diagonal= (elements+1)/2;

for i=1:diagonal
    fprintf(fileID,'%d, %d\n', i, i+1); %% node 1&2 of element
end

%truss members at the bottom
rest= elements-diagonal;
for i=1:rest
   fprintf(fileID, '%d, %d\n', i, i+2);
end

%truss members at the top
%top= bottom-1;
%for i=1:2:top
%   fprintf(fileID, '%d, %d\n', i+1, i+3);
%end

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



% for i=1:nodes
%    if i==1
% 	for j=1:3
% 	   fprintf(fileID, '%d, %d\n', i, j);
% 	end
%    end
%    if mod(i,2)==1
% 	fprintf(fileID, '%d, %d\n', i, 3);
%    end
%    if i==nodes
% 	for j=1:2:3
% 	  fprintf(fileID, '%d, %d\n', i, j );
% 	end
%    end
% end 
% 


%done with constraint statement:
fprintf(fileID,'0, 0\n');


xoffset=5;
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
fprintf(fileID,'2, 0.1, 0.1 200 0.001, 0.001');
fclose(fileID);
