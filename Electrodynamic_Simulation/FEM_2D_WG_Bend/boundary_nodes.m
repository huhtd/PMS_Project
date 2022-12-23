% Boundary nodes


% Vertical ports
xp1=0.;
yp2=0.12;

node_angle=zeros(1,Nn);

for i=1:Ne
    
    n1=el_no(i,1);
    n2=el_no(i,2);
    n3=el_no(i,3);
    
    x1=nodes(n1,1);
    y1=nodes(n1,2);
    x2=nodes(n2,1);
    y2=nodes(n2,2);
    x3=nodes(n3,1);
    y3=nodes(n3,2);
    
    angle1=acos(((x2-x1)*(x3-x1)+(y2-y1)*(y3-y1))/(sqrt((x2-x1)^2+(y2-y1)^2)...
        *sqrt((x3-x1)^2+(y3-y1)^2)));
    angle2=acos(((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2))/(sqrt((x1-x2)^2+(y1-y2)^2)...
        *sqrt((x3-x2)^2+(y3-y2)^2)));
    angle3=acos(((x1-x3)*(x2-x3)+(y1-y3)*(y2-y3))/(sqrt((x1-x3)^2+(y1-y3)^2)...
        *sqrt((x2-x3)^2+(y2-y3)^2)));

    node_angle(n1)=node_angle(n1)+abs(angle1);
    node_angle(n2)=node_angle(n2)+abs(angle2);
    node_angle(n3)=node_angle(n3)+abs(angle3);
    
end % i

for i=1:Nb
    
    n1=bc_elements(i,1);
    n2=bc_elements(i,2);
    
    x1=nodes(n1,1);
    y1=nodes(n1,2);
    x2=nodes(n2,1);
    y2=nodes(n2,2);
    
    if(abs(x1-x2)<0.0001 && abs(x1-xp1)<0.0001) % Input
        st_no(n1)=0;
        st_no(n2)=0;
    end
    if(abs(y1-y2)<0.0001 && abs(y1-yp2)<0.0001) % Output
        if (abs(x1-0.11)>0.0001 && abs(x1-0.13)>0.0001)
            st_no(n1)=0;
        end
        if (abs(x2-0.11)>0.0001 && abs(x2-0.13)>0.0001)
            st_no(n2)=0;
        end
    end
    
    if (abs(node_angle(n1)-2*pi))<0.0001 % not boundary node
            st_no(n1)=0;
    end
    if (abs(node_angle(n2)-2*pi))<0.0001 % not boundary node
            st_no(n2)=0;
    end
    
    if abs(y1-d/2)<0.0001 || abs(y1+d/2)<0.0001
        st_no(n1)=1;
    end
    
    if abs(y2-d/2)<0.0001 || abs(y2+d/2)<0.0001
        st_no(n2)=1;
    end
    
    
end % i
