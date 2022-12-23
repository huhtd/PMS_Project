% Poynting field computation
clear Sx;
clear Sy;
clear Smod;

Sx=zeros(1,Nn);
Sy=zeros(1,Nn);
Smod=zeros(1,Nn);

for i=1:Nn
   Sx(i)=-0.5*Ez(i)*conj(Hy(i)); 
   Sy(i)=0.5*Ez(i)*conj(Hx(i)); 
   Smod(i)=sqrt(Sx(i)^2+Sy(i)^2);
end

