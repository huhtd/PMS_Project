% Load and prepare data
clc;
clear;

load nodes.txt;
load elements.txt;
load domains.txt;
load bcs.txt;

[Nn m]=size(nodes);
[Ne m]=size(elements);
[Nb m]=size(bcs);

for i=1:Ne
   for j=1:3
       el_no(i,j)=elements(i,j)+1;  
   end
   el_mat(i)=domains(i);
end

for i=1:Nn
   x_no(i)=nodes(i,1);
   y_no(i)=nodes(i,2);
   st_no(i)=0;
end

for i=1:Nb
   bc_elements(i,1)=bcs(i,1)+1;
   bc_elements(i,2)=bcs(i,2)+1;
   in=bcs(i,1)+1;
   st_no(in)=1;
   in=bcs(i,2)+1;
   st_no(in)=1;
end

