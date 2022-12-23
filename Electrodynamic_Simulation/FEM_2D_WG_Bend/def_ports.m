% Assign BCS
% Ports definition

% verical ports
xp1=0.; % Input
yp2=0.12; % Output

imj=sqrt(-1);
for i=1:Nb % Loop over boundary line elements
    in1=bc_elements(i,1);
    in2=bc_elements(i,2);
    x1=x_no(in1);x2=x_no(in2);
    y1=y_no(in1);y2=y_no(in2);
    ysr=(y1+y2)/2;
    length=sqrt((x1-x2)^2+(y1-y2)^2);
    if(abs(x1-x2)<0.0001 && abs(x1-xp1)<0.0001) % Input
        A(in1,in2)=A(in1,in2)-imj*kx*length/6.;
        A(in2,in1)=A(in2,in1)-imj*kx*length/6.;
        A(in1,in1)=A(in1,in1)-imj*kx*length/3.;
        A(in2,in2)=A(in2,in2)-imj*kx*length/3.;
        b(in1)=b(in1)-2*imj*kx*length/2.*te_even_mode(ysr,E0,Order,d);
        b(in2)=b(in2)-2*imj*kx*length/2.*te_even_mode(ysr,E0,Order,d);
    end
    if(abs(y1-y2)<0.0001 && abs(y1-yp2)<0.0001) % Output
        A(in1,in2)=A(in1,in2)-imj*kx*length/6.;
        A(in2,in1)=A(in2,in1)-imj*kx*length/6.;
        A(in1,in1)=A(in1,in1)-imj*kx*length/3.;
        A(in2,in2)=A(in2,in2)-imj*kx*length/3.;
    end
end

return;
