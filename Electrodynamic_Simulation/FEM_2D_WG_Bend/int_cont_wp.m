% integral contribution of element
% function [ke]=int_cont_es(ie)
% ie - index of elemnt

in1=el_no(ie,1);
in2=el_no(ie,2);
in3=el_no(ie,3);

x1=x_no(in1);
x2=x_no(in2);
x3=x_no(in3);

y1=y_no(in1);
y2=y_no(in2);
y3=y_no(in3);

Se=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2.;
if(Se<0) 
    Se
end

if(el_mat(ie)==1)
    eps=eps1;
end
if(el_mat(ie)==2)
    eps=eps2;
end
if(el_mat(ie)==3)
    eps=eps3;
end

ae(1,1)=(y3-y2)^2+(x3-x2)^2;
ae(2,2)=(y3-y1)^2+(x3-x1)^2;
ae(3,3)=(y2-y1)^2+(x2-x1)^2;

ae(1,2)=(y3-y2)*(y1-y3)+(x3-x2)*(x1-x3);
ae(2,1)=ae(1,2);
ae(1,3)=(y2-y1)*(y3-y2)+(x2-x1)*(x3-x2);
ae(3,1)=ae(1,3);
ae(2,3)=(y1-y3)*(y2-y1)+(x1-x3)*(x2-x1);
ae(3,2)=ae(2,3);

for i=1:3
    for j=1:3
        ae(i,j)=-ae(i,j)/(4*Se);
        if (i==j)
            ae(i,j)=ae(i,j)+eps*k0^2*Se/6.;
        else
            ae(i,j)=ae(i,j)+eps*k0^2*Se/12.;
        end
    end
end
