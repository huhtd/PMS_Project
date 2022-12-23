% S11 - input port

% verical ports
xp1=0.; % Input
% xp2=0.22; % Output

S11=0;

integral1=0;
integral2=0;

imj=sqrt(-1);
for i=1:Nb % Loop over boundary line elements
    
    in1=bc_elements(i,1);
    in2=bc_elements(i,2);
    x1=x_no(in1);x2=x_no(in2);
    y1=y_no(in1);y2=y_no(in2);
    ysr=(y1+y2)/2;
    length=sqrt((x1-x2)^2+(y1-y2)^2);
    CF=0.5*(Ez(in1)+Ez(in2));
    MF=te_even_mode(ysr,E0,Order,d);
    if(abs(x1-x2)<0.0001 && abs(x1-xp1)<0.0001) % Input
        integral1=integral1+(CF-MF)*MF*length;
        integral2=integral2+MF*MF*length;
    end

end

S11=integral1/integral2;

return;
