% electric field computation and visualisation

% clear xe;
% clear ye;
clear Hx;
clear Hy;
clear Habs;
clear Ssupp;

Hx=zeros(1,Nn);
Hy=zeros(1,Nn);
Habs=zeros(1,Nn);
Ssupp=zeros(1,Nn);

ind=0;

for i=1:Ne
    %if(el_mat(i)==1)
    %    eps=eps1;
    %else
    %    eps=eps2;
    %end
    
    in1=el_no(i,1);
    in2=el_no(i,2);
    in3=el_no(i,3);
    
    V1=Vc(in1);
    V2=Vc(in2);
    V3=Vc(in3);
    
    x1=x_no(in1);
    y1=y_no(in1);
    x2=x_no(in2);
    y2=y_no(in2);
    x3=x_no(in3);
    y3=y_no(in3);
    
    Se=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2.;
        
    %Elemental values
    Hey=1/(imj*omega*mu0)*((y2-y3)*V1+(y3-y1)*V2+(y1-y2)*V3)/(2*Se);
    Hex=-1/(imj*omega*mu0)*((x3-x2)*V1+(x1-x3)*V2+(x2-x1)*V3)/(2*Se);
    
    %Averaging to nodes    
    Hy(in1)=Hy(in1)+Se*(Hey);
    Hx(in1)=Hx(in1)+Se*(Hex);
    Ssupp(in1)=Ssupp(in1)+Se;
    
    Hy(in2)=Hy(in2)+Se*(Hey);
    Hx(in2)=Hx(in2)+Se*(Hex);
    Ssupp(in2)=Ssupp(in2)+Se;
    
    Hy(in3)=Hy(in3)+Se*(Hey);
    Hx(in3)=Hx(in3)+Se*(Hex);
    Ssupp(in3)=Ssupp(in3)+Se;
    
end

for i=1:Nn
    Hx(i)=Hx(i)/Ssupp(i);
    Hy(i)=Hy(i)/Ssupp(i);
    Habs(i)=sqrt(Hx(i)^2+Hy(i)^2);
end

return;

figure(3);
x_min=min(x_no);
x_max=max(x_no);
y_min=min(y_no);
y_max=max(y_no);
axis([x_min x_max y_min y_max]);
hold on;
axis equal;
axis off;

quiver(x_no,y_no,real(Hx),real(Hy),1,'k');


